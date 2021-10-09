#include "DifferentiallyPrivateNonparametricHypthothesisTesting/privatekruskalwallis/privatekruskalwallis.h"

double dpnht::internal_kruskal_wallis::KruskalWallisAbs(const std::vector<std::vector<double>>& rankings_data) {
	// compute average of all ranks
	double average_rank = 0;
	double num_groups = rankings_data.size();
	for (int i = 0; i < num_groups; i++) {
		average_rank += rankings_data[i].size();
	}
	average_rank = (average_rank + 1) / 2;

	double numerator = 0;
	double denominator = 0;

	// compute numerator and denominator simultaneously
	for (int i = 0; i < num_groups; i++) {
		double running_total = 0;
		int group_size = rankings_data[i].size();
		for (int j = 0; j < group_size; j++) {
			double temp = rankings_data[i][j];
			running_total += temp;
			denominator += std::abs(temp - average_rank);
		}
		numerator += group_size * std::abs((running_total / group_size) - average_rank);
	}

	return (num_groups - 1) * numerator - denominator;
}