#include "miner.h"
#include "utils.h"
#include "neighborhood_mgr.h"
#include "types.h"
#include <algorithm>
#include <set>
#include <string>
#include <iostream>
#include <omp.h> 
#include <iomanip>
#include <chrono>

// -----------------------------------------------------------------------------
// Main Mining Function
// -----------------------------------------------------------------------------
std::vector<std::pair<Colocation, double>> MaxPRMiner::mineColocations(
    double minMaxPR,
    NeighborhoodMgr* neighborhoodMgr,
    const std::vector<SpatialInstance>& instances
) {
	// 1. Preprocessing: count total instances for each feature type (global count)
    globalFeatureCounts.clear();
    for (const auto& inst : instances) {
        globalFeatureCounts[inst.type]++;
    }

    // Result container
    std::vector<std::pair<Colocation, double>> resultPatterns;

    // Store Pattern and Rowset for the current k step
    // Key: Pattern (sorted), Value: {maxPR, vector<RowInstance>}
    std::map<Colocation, PatternData> prevMap;

    // -------------------------------------------------------
    // Step k = 2: Generate from Neighborhood Graph (Geometric Method)
    // -------------------------------------------------------
    int k = 2;
    // Temporary map to group instance pairs into the correct pattern
    std::map<Colocation, std::vector<ColocationInstance>> c2Candidates;

    auto allStars = neighborhoodMgr->getAllStarNeighborhoods();
    for (const auto& typePair : allStars) {
        for (const auto& star : typePair.second) {
            const SpatialInstance* centerInst = star.center;
            for (const auto* neighborInst : star.neighbors) {
                if (centerInst->type >= neighborInst->type) continue;
                Colocation pattern = { centerInst->type, neighborInst->type };
                ColocationInstance row = { centerInst, neighborInst };

                c2Candidates[pattern].push_back(row);
            }
        }
    }

	// filter candidates k=2 by MaxPR and prepare for next step
    std::vector<Colocation> prevPatternsList;
    std::set<Colocation> prevPatternSet;

    for (auto& item : c2Candidates) {
        double pr = calculateMaxPR(item.first, item.second);
        if (pr >= minMaxPR) {
            prevMap[item.first] = { pr, std::move(item.second) };
            prevPatternsList.push_back(item.first);
            prevPatternSet.insert(item.first);
            resultPatterns.push_back({ item.first, pr });
        }
    }

    // -------------------------------------------------------
    // Step k > 2: Loop
    // -------------------------------------------------------
    while (!prevPatternsList.empty()) {
        k++;
        std::map<Colocation, PatternData> currentMap;
        std::vector<Colocation> candidates = generateCandidatesMaxPR(prevPatternsList, prevPatternSet, minMaxPR);

        if (candidates.empty()) break;
        for (const auto& candPattern : candidates) {
            // Split candidate k into 2 pattern size (k-1) to join
            // Assume candidate is already sorted: [0, 1, ..., k-1]
            // Sub1: [0, 1, ..., k-2] (remove last element)
            // Sub2: [0, 1, ..., k-3, k-1] (remove second last element)

            Colocation sub1 = candPattern;
            sub1.pop_back(); // Remove last element

            Colocation sub2 = candPattern;
            sub2.erase(sub2.begin() + (sub2.size() - 2)); // Remove second last element

            // Check if the 2 sub-patterns exist in prevMap
            // (Due to weak monotonicity, one of the sub-patterns may not exist, but here
            // the candidate generation logic ensures we join from 2 existing ones)

            auto it1 = prevMap.find(sub1);
            auto it2 = prevMap.find(sub2);

            if (it1 != prevMap.end() && it2 != prevMap.end()) {
                std::vector<ColocationInstance> newRowset = generateRowset(it1->second.rowset, it2->second.rowset, k, neighborhoodMgr);
                double pr = calculateMaxPR(candPattern, newRowset);

                if (pr >= minMaxPR) {
                    currentMap[candPattern] = { pr, std::move(newRowset) };
                    resultPatterns.push_back({ candPattern, pr });
                }
            }
        }

		// 3. Prepare for next iteration
		prevMap = std::move(currentMap); //Delete old map to free memory

        prevPatternsList.clear();
        prevPatternSet.clear();
        for (const auto& item : prevMap) {
            prevPatternsList.push_back(item.first);
            prevPatternSet.insert(item.first);
        }
    }

    return resultPatterns;
}

// -----------------------------------------------------------------------------
// Helper: Calculate MaxPR for a given pattern and its rowset
// Formula: max( PR(f1), PR(f2)... )
// PR(f) = (Number of unique instances of f participating in the pattern) / (Total instances of f)
// -----------------------------------------------------------------------------
double MaxPRMiner::calculateMaxPR(const Colocation& pattern, const std::vector<ColocationInstance>& rowset) {
    if (rowset.empty()) return 0.0;

    double maxVal = 0.0;
    for (size_t i = 0; i < pattern.size(); ++i) {
        FeatureType type = pattern[i];
        std::set<const SpatialInstance*> uniqueInstances;
        for (const auto& row : rowset) {
            uniqueInstances.insert(row[i]);
        }

        int countInPattern = uniqueInstances.size();
        int totalGlobal = globalFeatureCounts[type];

        if (totalGlobal > 0) {
            double pr = (double)countInPattern / totalGlobal;
            if (pr > maxVal) maxVal = pr;
        }
    }
    return maxVal;
}

// -----------------------------------------------------------------------------
// Helper: Generate Rowset (Join Step)
// Concat rows from rows1 and rows2 based on the join condition (k-2 prefix match) and clique condition (new edge)
// Assumption: rows1 has the form {x1...x(k-2), A}, rows2 has the form {x1...x(k-2), B}
// Result: {x1...x(k-2), A, B} if A and B are neighbors
// -----------------------------------------------------------------------------
std::vector<ColocationInstance> MaxPRMiner::generateRowset(
    const std::vector<ColocationInstance>& rows1,
    const std::vector<ColocationInstance>& rows2,
    int k,
    NeighborhoodMgr* nbrMgr
) {
    std::vector<ColocationInstance> result;
    result.reserve(std::min(rows1.size(), rows2.size()));

    // Lưu ý: rows1 và rows2 thường đã được sort theo các phần tử đầu tiên (do cách sinh candidate)
    // Nếu chưa sort, cần sort hoặc dùng hash map để join nhanh.
    // Ở đây dùng giải thuật Nested Loop Join đơn giản (có thể tối ưu thành Sort-Merge Join).

    // Tối ưu: Xây dựng Map cho rows2 dựa trên (k-2) phần tử đầu
    // Key: Vector<Instance*> (prefix), Value: Vector<Instance*> (phần tử cuối)
    // Tuy nhiên, để đơn giản và hiệu quả bộ nhớ, ta dùng 2 con trỏ nếu dữ liệu đã sort, 
    // hoặc Multimap cho prefix.

    // Dưới đây là cách cài đặt Sort-Merge Join giả định các rowset đã sort theo các instance ID
    // Để an toàn, ta dùng Multimap cho rows2 để lookup nhanh

    // Key: ID của (k-2) instance đầu tiên. Value: index trong rows2
    // Để key đơn giản, dùng ID instance đầu tiên (hoặc hash của prefix)
    using Prefix = std::vector<const SpatialInstance*>;
    std::map<Prefix, std::vector<size_t>> mapRows2;

    for (size_t i = 0; i < rows2.size(); ++i) {
        Prefix prefix(rows2[i].begin(), rows2[i].end() - 1);
        mapRows2[prefix].push_back(i);
    }

    for (const auto& r1 : rows1) {
        Prefix prefix1(r1.begin(), r1.end() - 1);

        auto it = mapRows2.find(prefix1);
        if (it != mapRows2.end()) {
            const SpatialInstance* instA = r1.back();

            for (size_t idx2 : it->second) {
                const SpatialInstance* instB = rows2[idx2].back();
                if (instA == instB) continue;
                if (nbrMgr->areNeighbors(instA, instB)) {
                    ColocationInstance newRow = r1;
                    newRow.push_back(instB);
                    result.push_back(newRow);
                }
            }
        }
    }

    return result;
}

// -----------------------------------------------------------------------------
// Logic generate candidate patterns based on Weak Monotonicity
// -----------------------------------------------------------------------------
std::vector<Colocation> MaxPRMiner::generateCandidatesMaxPR(
    const std::vector<Colocation>& prevMaxPR,
    const std::set<Colocation>& prevSet,
    double minMaxPR)
{
    std::vector<Colocation> candidates;
    if (prevMaxPR.empty()) return candidates;

    size_t k_minus_1 = prevMaxPR[0].size();

    // Self-join prevMaxPR
    for (size_t i = 0; i < prevMaxPR.size(); ++i) {
        for (size_t j = i + 1; j < prevMaxPR.size(); ++j) {
            const auto& p1 = prevMaxPR[i];
            const auto& p2 = prevMaxPR[j];

			// 1. Check join condition: p1 and p2 must share the same prefix of length k-2
            // p1: [a, b, X]
            // p2: [a, b, Y]
            bool match = true;
            for (size_t x = 0; x < k_minus_1 - 1; ++x) {
                if (p1[x] != p2[x]) {
                    match = false; break;
                }
            }
            if (!match) continue; 

			// 2. Create new candidate by merging p1 and p2
            Colocation candidate = p1;
            candidate.push_back(p2.back());
            if (p1.back() >= p2.back()) continue;

			// 3. Pruning using Weak Monotonicity: Check if candidate can be accepted based on its (k-1) subsets
            // A k-candidate is accepted if it has AT MOST 1 (k-1) subset not in prevSet (i.e., < threshold)

            int missingSubsets = 0;
            for (size_t skip = 0; skip < candidate.size(); ++skip) {
                Colocation subset = candidate;
                subset.erase(subset.begin() + skip);
                if (prevSet.find(subset) == prevSet.end()) {
                    missingSubsets++;
                }
                if (missingSubsets > 1) break;
            }

            if (missingSubsets <= 1) {
                candidates.push_back(candidate);
            }
        }
    }
    return candidates;
}