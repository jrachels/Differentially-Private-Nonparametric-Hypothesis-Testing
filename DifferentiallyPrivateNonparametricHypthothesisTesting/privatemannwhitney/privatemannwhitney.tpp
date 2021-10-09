#ifndef DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEMANNWHITNEY_TPP
#define DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEMANNWHITNEY_TPP

#include <iostream>

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/privatemannwhitney/privatemannwhitney.h"

template<std::uniform_random_bit_generator Engine>
std::pair<double, double> dpnht::internal_mann_whitney::PrivateMannWhitneyStatistic(const std::pair<std::vector<double>, std::vector<double>>& rankings_data, const double epsilon_m, const double epsilon_U, const double delta, Engine& Eng) {
	std::pair<int, int> size_of_groups(rankings_data.first.size(), rankings_data.second.size());
	double m = std::min(size_of_groups.first, size_of_groups.second);
	laplace_distribution  ld_m{ 0, 1 / epsilon_m };
	double m_approx = m + ld_m(Eng);

	double c = std::log(2 * delta) / epsilon_m;

	double m_star = std::max<double>(std::ceil(m_approx - c), 0);


	// POSSIBLE ERROR HERE: I AM FORCED TO ROUND UP TO KEEP SCALE ARGUMENT OF LAPLACE POSITIVE

	laplace_distribution  ld_U{ 0, std::max<double>((size_of_groups.first + size_of_groups.second - m_star), 0) / epsilon_U };

	double private_u_statistic = PublicMannWhitney(rankings_data) + ld_U(Eng);

	std::pair<double, double> mann_whitney_statistics(m_approx, private_u_statistic);

	return mann_whitney_statistics;
}

template<std::uniform_random_bit_generator Engine>
std::pair<std::vector<double>, std::vector<double>> dpnht::internal_mann_whitney::BuildRandomDatabase(const std::pair<int, int>& group_sizes, Engine& Eng) {

	// do some assertions

	std::pair<std::vector<double>, std::vector<double>> raw_data;

	raw_data.first.resize(group_sizes.first);
	raw_data.first.shrink_to_fit();
	for (std::size_t i = 0; i < raw_data.first.size(); i++) {
		raw_data.first[i] = std::generate_canonical<double, static_cast<size_t>(std::numeric_limits<double>::digits)>(Eng);
	}

	raw_data.second.resize(group_sizes.second);
	raw_data.second.shrink_to_fit();
	for (std::size_t i = 0; i < raw_data.second.size(); i++) {
		raw_data.second[i] = std::generate_canonical<double, static_cast<size_t>(std::numeric_limits<double>::digits)>(Eng);
	}

	return raw_data;
}


// TODO: Maybe rankdata shouldn't store the datatype unless the datatype takes a long time to copy
template<std::totally_ordered datatype>
std::pair<std::vector<double>, std::vector<double>> dpnht::internal_mann_whitney::ComputeRankings(const std::pair<std::vector<datatype>, std::vector<datatype>>& raw_data) {
	std::pair<std::vector<double>, std::vector<double>> rank_data;
	// 1. create map. keys are data point, values are vectors of columns containing that data point
	std::map<datatype, std::pair<int, int>> raw_data_tree;
	// 2. fill up map with raw_data

	for (std::size_t j = 0; j < raw_data.first.size(); j++) {
		auto [iter, success] = raw_data_tree.try_emplace(raw_data.first[j], 1, 0);
		if (!success) {
			(iter->second).first += 1;
		}
	}
	for (std::size_t j = 0; j < raw_data.second.size(); j++) {
		auto [iter, success] = raw_data_tree.try_emplace(raw_data.second[j], 0, 1);
		if (!success) {
			(iter->second).second += 1;
		}
	}
	// 3. iterate through map to fill up rankings data with rankings
	double rank = 0;

	for (auto iter = raw_data_tree.begin(); iter != raw_data_tree.end(); ++iter) {
		int num_first_group = (iter->second).first;
		int num_second_group = (iter->second).second;
		double num_ties = num_first_group + num_second_group; // conversion here is fine overhead.
		double average_rankings = rank + (num_ties + 1) / 2.0;
		for (int i = 0; i < num_first_group; i++) {
			rank_data.first.push_back(average_rankings);
		}
		for (int i = 0; i < num_second_group; i++) {
			rank_data.second.push_back(average_rankings);
		}
		rank += num_ties;
	}
	return rank_data;
}

template<std::uniform_random_bit_generator Engine>
double dpnht::internal_mann_whitney::ComputePValue(const std::pair<int, int>& group_sizes, const double private_u_statistic, const double epsilon_m, const double epsilon_u, const double delta, Engine& Eng, const int p_value_trials) {
	// assert some things to make sure division and mod works as expected
	int successes = 0;

	for (int i = 0; i < p_value_trials; i++) {
		std::pair<std::vector<double>, std::vector<double>> raw_data = BuildRandomDatabase(group_sizes, Eng);
		auto [m_approx_trial, private_u_static_trial] = PrivateMannWhitneyStatistic(ComputeRankings(raw_data), epsilon_m, epsilon_u, delta, Eng);
		if (private_u_statistic > private_u_static_trial) {
			successes += 1;
		}
	}
	return static_cast<double>(successes) / p_value_trials;
}

template<std::totally_ordered datatype, std::uniform_random_bit_generator Engine>
std::pair<const double, const double> dpnht::internal_mann_whitney::PrivateMannWhitneyHelper(const std::pair<std::vector<datatype>, std::vector<datatype>>& raw_data, const double epsilon_m, const double epsilon_u, const double delta, Engine& Eng, const int p_value_trials) {
	// compute ranks, handling ties
	std::pair<std::vector<double>, std::vector<double>> rankings_data = ComputeRankings(raw_data);
	int num_data_points = rankings_data.first.size() + rankings_data.second.size();
	// find statistic
	auto [m_approx, private_u_statistic] = PrivateMannWhitneyStatistic(rankings_data, epsilon_m, epsilon_u, delta, Eng);
	int m_approx_rounded = static_cast<int>(std::max<double>(0, m_approx));


	// POSSIBLE ERROR HERE: I AM FORCED TO CAP M_APPROX TO KEEP GROUP SIZES POSITIVE

	m_approx_rounded = std::min<int>(m_approx_rounded, num_data_points);
	// compute p value
	std::pair<int, int> group_sizes(m_approx_rounded, num_data_points - m_approx_rounded);
	double p_value = ComputePValue(group_sizes, private_u_statistic, epsilon_m, epsilon_u, delta, Eng, p_value_trials);

	// return statistic and p value
	return std::pair<const double, const double>(private_u_statistic, p_value);
}

template<std::totally_ordered datatype, std::uniform_random_bit_generator Engine>
[[nodiscard]] std::pair<const double, const double> dpnht::PrivateMannWhitney(const std::pair<std::vector<datatype>, std::vector<datatype>>& raw_data, const double epsilon_m, const double epsilon_u, const double delta, Engine& Eng, const int p_value_trials) {
	return dpnht::internal_mann_whitney::PrivateMannWhitneyHelper(raw_data, epsilon_m, epsilon_u, delta, Eng, p_value_trials);
}

template<std::totally_ordered datatype>
[[nodiscard]] std::pair<const double, const double> dpnht::PrivateMannWhitney(const std::pair<std::vector<datatype>, std::vector<datatype>>& raw_data, const double epsilon_m, const double epsilon_u, const double delta, const int p_value_trials) {
	// create random engine
	std::cout << "Warning: PrivateMannWhitney is being called without a random engine. An std::mt19937 will be instantiated with random_device() instead. See comments for details.\n";
	std::mt19937 mt{ std::random_device()() };
	// Exceptions
	return dpnht::internal_mann_whitney::PrivateMannWhitneyHelper(raw_data, epsilon_m, epsilon_u, delta, mt, p_value_trials);

}

#endif