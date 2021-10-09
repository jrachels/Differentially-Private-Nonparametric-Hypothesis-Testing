#include "DifferentiallyPrivateNonparametricHypthothesisTesting/privatemannwhitney/privatemannwhitney.h"

double dpnht::internal_mann_whitney::PublicMannWhitney(const std::pair<std::vector<double>, std::vector<double>>& rankings_data) {
	// compute average of all ranks
	std::pair<double, double> u_statistic_intermediates(0, 0);
	std::pair<int, int> size_of_groups(rankings_data.first.size(), rankings_data.second.size());
	for (int i = 0; i < size_of_groups.first; i++) {
		u_statistic_intermediates.first += rankings_data.first[i];
	}
	for (int j = 0; j < size_of_groups.second; j++) {
		u_statistic_intermediates.second += rankings_data.second[j];
	}
	u_statistic_intermediates.first -= size_of_groups.first * (size_of_groups.first + 1) / 2;
	u_statistic_intermediates.second -= size_of_groups.second * (size_of_groups.second + 1) / 2;

	return std::min(u_statistic_intermediates.first, u_statistic_intermediates.second);
}