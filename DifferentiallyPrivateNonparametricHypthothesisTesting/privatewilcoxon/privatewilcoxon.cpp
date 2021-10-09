#include <iostream>

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/privatewilcoxon/privatewilcoxon.h"

// Pratt Variant
double dpnht::internal_wilcoxon::PublicWilcoxon(const std::vector<std::pair<double, double>>& dataset) {
	double public_w_statistic = 0;
	// key is absolute difference. value is (sum  of signs, number of datapoints) pairs.
	std::map<double, std::pair<int, int>> data_set_tree;


	int size = dataset.size();


	for (int i = 0; i < size; i++) {
		double u_i = dataset[i].first;
		double v_i = dataset[i].second;

		double absolute_difference = v_i - u_i; // may not be positive yet
		int sign = 0;
		if (absolute_difference > 0) {
			sign = 1;
		}
		else if (absolute_difference < 0) {
			absolute_difference = -absolute_difference;
			sign = -1;
		}

		auto [iter, success] = data_set_tree.try_emplace(absolute_difference, sign, 1);
		if (!success) {
			(iter->second).first += sign;
			(iter->second).second += 1;
		}
	}

	// process public_w_statistic
	double rank = 0;

	for (auto iter = data_set_tree.begin(); iter != data_set_tree.end(); ++iter) {
		double num_ties = (iter->second).second;
		double average_rankings = rank + (num_ties + 1) / 2.0;
		public_w_statistic += average_rankings * ((iter->second).first);
		rank += num_ties;
	}

	assert(rank == size);
	return public_w_statistic;

}


[[nodiscard]] std::pair<const double, const double> dpnht::PrivateWilcoxon(const std::vector<std::pair<double, double>>& dataset, const double epsilon, const int p_value_trials) {
	// create random engine
	std::cout << "Warning: PrivateWilcoxon is being called without a random engine. An std::mt19937 will be instantiated with random_device() instead. See comments for details.\n";
	std::mt19937 mt{ std::random_device()() };
	// Exceptions
	return dpnht::internal_wilcoxon::PrivateWilcoxonHelper(dataset, epsilon, mt, p_value_trials);

}