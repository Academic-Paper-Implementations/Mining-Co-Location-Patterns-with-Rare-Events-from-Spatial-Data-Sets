#pragma once
#include "types.h"
#include "neighborhood_mgr.h"
#include <vector>
#include <map>
#include <functional>
#include <set>


struct PatternData {
    double maxPR;
    std::vector<ColocationInstance> rowset;
};

class MaxPRMiner {
private:
	// count number of instances for each feature type in the entire dataset (global count)
    std::map<FeatureType, int> globalFeatureCounts;

	// Helper: Calculate MaxPR for a given pattern and its rowset
    double calculateMaxPR(const Colocation& pattern, const std::vector<ColocationInstance>& rowset);

    // Helper: Generate Rowset for candidate k from 2 parent (k-1)
    std::vector<ColocationInstance> generateRowset(
        const std::vector<ColocationInstance>& rows1,
        const std::vector<ColocationInstance>& rows2,
        int k,
        NeighborhoodMgr* nbrMgr
    );

public:
    std::vector<std::pair<Colocation, double>> mineColocations(
        double minMaxPR,
        NeighborhoodMgr* nbrMgr,
        const std::vector<SpatialInstance>& instances
    );
    // Logic to generate candidate patterns based on Weak Monotonicity
    std::vector<Colocation> generateCandidatesMaxPR(
        const std::vector<Colocation>& prevPatterns,
        const std::set<Colocation>& prevPatternSet, // Pass Set for O(logN) lookup
        double minMaxPR
    );
};