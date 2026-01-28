/**
 * @file spatial_index.cpp
 * @brief Implementation of spatial indexing and neighbor search
 */

#include "spatial_index.h"
#include <cmath>
#include <iostream>
#include <algorithm>

SpatialIndex::SpatialIndex(double distThresh)
    : distanceThreshold(distThresh)
{

}

double SpatialIndex::euclideanDist(const SpatialInstance& a, const SpatialInstance& b) {
    // Calculate Euclidean distance using the Pythagorean theorem
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

std::vector<std::pair<SpatialInstance, SpatialInstance>> SpatialIndex::findNeighborPair(const std::vector<SpatialInstance>& instances) {
    std::vector<std::pair<SpatialInstance, SpatialInstance>> neighborPairs;

	//Divide to grid cells for optimization
    double minX = std::min_element(instances.begin(), instances.end(),
		[](const SpatialInstance& a, const SpatialInstance& b) { return a.x < b.x; })->x;
	double minY = std::min_element(instances.begin(), instances.end(),
		[](const SpatialInstance& a, const SpatialInstance& b) { return a.y < b.y; })->y;
	double maxX = std::max_element(instances.begin(), instances.end(),
		[](const SpatialInstance& a, const SpatialInstance& b) { return a.x < b.x; })->x;
	double maxY = std::max_element(instances.begin(), instances.end(),
		[](const SpatialInstance& a, const SpatialInstance& b) { return a.y < b.y; })->y;

	int gridX = static_cast<int>(std::ceil((maxX - minX) / distanceThreshold));
	int gridY = static_cast<int>(std::ceil((maxY - minY) / distanceThreshold));

	std::vector<std::vector<SpatialInstance>> gridCells(gridX * gridY);

    for (const auto& inst : instances) {
        int cx = static_cast<int>((inst.x - minX) / distanceThreshold);
        int cy = static_cast<int>((inst.y - minY) / distanceThreshold);

        gridCells[cx * gridY + cy].push_back(inst);
    }

    for (int cx = 0; cx < gridX; ++cx) {
        for (int cy = 0; cy < gridY; ++cy) {

            const auto& cell = gridCells[cx * gridY + cy];

            // compare between instance in center cell with around cells
            for (size_t i = 0; i < cell.size(); ++i) {
				for (size_t j = i + 1; j < cell.size(); ++j) {
                    if (cell[i].type != cell[j].type && euclideanDist(cell[i], cell[j]) <= distanceThreshold) {
                        neighborPairs.emplace_back(cell[i], cell[j]);
                    }
                }
                for (int dx = 0; dx <= 1; ++dx) {
                    for (int dy = (dx == 0 ? 1 : -1); dy <= 1; ++dy) {
                        if (dx == 0 && dy == 0) continue; // skip center cell
                        int nx = cx + dx;
                        int ny = cy + dy;
                        if (nx >= 0 && nx < gridX && ny >= 0 && ny < gridY) {
                            const auto& neighborCell = gridCells[nx * gridY + ny];
                            for (const auto& neighborInst : neighborCell) {
                                if (cell[i].type != neighborInst.type && euclideanDist(cell[i], neighborInst) <= distanceThreshold) {
                                    neighborPairs.emplace_back(cell[i], neighborInst);
                                }
                            }
                        }
                    }
				}
            }
        }
    }
    for (auto& p : neighborPairs) {
        bool needSwap = false;

        // So sánh Feature (String) trước
        if (p.first.type > p.second.type) {
            needSwap = true;
        }
        // Nếu Feature giống nhau thì so sánh ID
        else if (p.first.type == p.second.type && p.first.id > p.second.id) {
            needSwap = true;
        }

        if (needSwap) {
            std::swap(p.first, p.second);
        }
    }

    return neighborPairs;
}