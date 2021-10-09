#ifndef DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEKRUSKALWALLIS_TPP
#define DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEKRUSKALWALLIS_TPP

#include <map>
#include <iostream>

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/privatekruskalwallis/privatekruskalwallis.h"

template<std::uniform_random_bit_generator Engine>
double dpnht::internal_kruskal_wallis::PrivateKruskalWallisStatistic<Engine>(const std::vector<std::vector<double>>& rankings_data, double epsilon, Engine& Eng) {
	laplace_distribution  ld{ 0, 8 / epsilon };
	double public_h_statistic = KruskalWallisAbs(rankings_data);
	return public_h_statistic + ld(Eng);
}



template<std::uniform_random_bit_generator Engine>
std::vector<std::vector<double>> dpnht::internal_kruskal_wallis::BuildRandomDatabase<Engine>(const std::size_t num_groups, const std::size_t num_data_points, Engine& Eng) {

	// do some assertions

	std::vector<std::vector<double>> raw_data(num_groups);
	const std::size_t size_groups = num_data_points / num_groups;
	const std::size_t left_over = num_data_points % num_groups;
	for (std::size_t i = 0; i < left_over; i++) {
		raw_data[i].resize(size_groups + 1);
		raw_data.shrink_to_fit();
		// fill this entry with random numbers
		for (std::size_t j = 0; j < raw_data[i].size(); j++) {
			raw_data[i][j] = std::generate_canonical<double, static_cast<size_t>(-1)>(Eng);
		}
	}
	for (std::size_t i = left_over; i < num_groups; i++) {
		raw_data[i].resize(size_groups);
		raw_data.shrink_to_fit();
		// fill this entry with random numbers
		for (std::size_t j = 0; j < raw_data[i].size(); j++) {
			raw_data[i][j] = std::generate_canonical<double, static_cast<size_t>(std::numeric_limits<double>::digits)>(Eng);
		}
	}
	return raw_data;
}

template<std::uniform_random_bit_generator Engine>
std::vector<std::vector<double>> dpnht::internal_kruskal_wallis::BuildRandomDatabase<Engine>(const std::vector<int>& group_sizes, Engine& Eng) {

	// do some assertions

	std::vector<std::vector<double>> raw_data(group_sizes.size());

	for (std::size_t i = 0; i < raw_data.size(); i++) {
		raw_data[i].resize(group_sizes[i]);
		raw_data[i].shrink_to_fit();
		for (std::size_t j = 0; j < raw_data[i].size(); j++) {
			raw_data[i][j] = std::generate_canonical<double, static_cast<size_t>(std::numeric_limits<double>::digits)>(Eng);
		}
	}

	return raw_data;
}

template<std::totally_ordered datatype>
std::vector<std::vector<double>> dpnht::internal_kruskal_wallis::ComputeRankings<datatype>(const std::vector<std::vector<datatype>>& raw_data) {
	std::vector<std::vector<double>> rank_data(raw_data.size());
	// 1. create map. keys are data point, values are vectors of columns containing that data point
	std::map<datatype, std::vector<std::size_t>> raw_data_tree;
	// 2. fill up map with raw_data
	for (std::size_t i = 0; i < raw_data.size(); i++) {
		for (std::size_t j = 0; j < raw_data[i].size(); j++) {
			auto [iter, success] = raw_data_tree.try_emplace(raw_data[i][j], std::initializer_list<std::size_t>{i});
			if (!success) {
				(iter->second).push_back(i);
			}
		}
	}
	// 3. iterate through map to fill up rankings data with rankings
	double rank = 0;

	//for (std::map<datatype, std::vector<int>>::iterator iter = raw_data_tree.begin(); iter != raw_data_tree.end(); ++iter) {

	// replace auto with the actual type

	for (auto iter = raw_data_tree.begin(); iter != raw_data_tree.end(); ++iter) {
		double num_ties = (iter->second).size();
		double average_rankings = rank + (num_ties + 1) / 2.0;
		for (int x : (iter->second)) {
			rank_data[x].push_back(average_rankings);
		}
		rank += num_ties;
	}

	return rank_data;
}

template<std::uniform_random_bit_generator Engine>
double dpnht::internal_kruskal_wallis::ComputePValue<Engine>(const std::vector<int>& group_sizes, const double private_h_statistic, const double epsilon, Engine& Eng, const int p_value_trials) {
	// assert some things to make sure division and mod works as expected
	int successes = 0;

	for (int i = 0; i < p_value_trials; i++) {
		std::vector<std::vector<double>> raw_data = BuildRandomDatabase(group_sizes, Eng);
		const double private_h_static_trial = PrivateKruskalWallisStatistic(ComputeRankings(raw_data), epsilon, Eng);
		if (private_h_statistic < private_h_static_trial) {
			successes += 1;
		}
	}
	return static_cast<double>(successes) / p_value_trials;
}

// Note: this provides a conservative estimate of the p-value without revealing group sizes
template<std::uniform_random_bit_generator Engine>
double dpnht::internal_kruskal_wallis::ComputePValue<Engine>(const int num_groups, const int num_data_points, const double private_h_statistic, const double epsilon, Engine& Eng, const int p_value_trials) {
	// assert some things to make sure division and mod works as expected
	int successes = 0;

	for (int i = 0; i < p_value_trials; i++) {
		std::vector<std::vector<double>> raw_data = BuildRandomDatabase(num_groups, num_data_points, Eng);
		const double private_h_static_trial = PrivateKruskalWallisStatistic(ComputeRankings(raw_data), epsilon, Eng);
		if (private_h_statistic < private_h_static_trial) {
			successes += 1;
		}
	}
	return static_cast<double>(successes) / p_value_trials;
}

template<std::totally_ordered datatype, std::uniform_random_bit_generator Engine>
[[nodiscard]] std::pair<const double, const double> dpnht::PrivateKruskalWallis<datatype, Engine>(const std::vector<std::vector<datatype>>& raw_data, const double epsilon, Engine& Eng, const int p_value_trials, bool group_sizes_known) {

	// Exceptions

	// compute ranks, handling ties
	std::vector<std::vector<double>> rankings_data = dpnht::internal_kruskal_wallis::ComputeRankings(raw_data);
	// find statistic
	double private_h_statistic = dpnht::internal_kruskal_wallis::PrivateKruskalWallisStatistic(rankings_data, epsilon, Eng);
	// compute p value
	double p_value = 0;
	if (group_sizes_known) {
		int num_groups = raw_data.size();
		std::vector<int> group_sizes(num_groups);
		group_sizes.shrink_to_fit();
		for (int i = 0; i < num_groups; i++) {
			group_sizes[i] = raw_data[i].size();
		}
		double p_value = dpnht::internal_kruskal_wallis::ComputePValue(group_sizes, private_h_statistic, epsilon, Eng, p_value_trials);
	}
	else {
		int num_groups = raw_data.size();
		int num_data_points = 0;
		for (int i = 0; i < num_groups; i++) {
			num_data_points += raw_data[i].size();
		}
		p_value = dpnht::internal_kruskal_wallis::ComputePValue(num_groups, num_data_points, private_h_statistic, epsilon, Eng, p_value_trials);
	}
	// return statistic and p value

	return std::pair<const double, const double>(private_h_statistic, p_value);
}

template<std::totally_ordered datatype>
[[nodiscard]] std::pair<const double, const double> dpnht::PrivateKruskalWallis(const std::vector<std::vector<datatype>>& raw_data, const double epsilon, const int p_value_trials, bool group_sizes_known) {
	// create random engine
	std::cout << "Warning: PrivateKruskalWallis is being called without a random engine. An std::mt19937 will be instantiated with random_device() instead. See comments for details.\n";
	std::mt19937 mt{ std::random_device()() };
	// Exceptions


	// compute ranks, handling ties
	std::vector<std::vector<double>> rankings_data = dpnht::internal_kruskal_wallis::ComputeRankings(raw_data);
	// find statistic
	double private_h_statistic = dpnht::internal_kruskal_wallis::PrivateKruskalWallisStatistic(rankings_data, epsilon, mt);
	// compute p value
	double p_value = 0;
	if (group_sizes_known) {
		int num_groups = raw_data.size();
		std::vector<int> group_sizes(num_groups);
		group_sizes.shrink_to_fit();
		for (int i = 0; i < num_groups; i++) {
			group_sizes[i] = raw_data[i].size();
		}
		double p_value = dpnht::internal_kruskal_wallis::ComputePValue(group_sizes, private_h_statistic, epsilon, mt, p_value_trials);
	}
	else {
		int num_groups = raw_data.size();
		int num_data_points = 0;
		for (int i = 0; i < num_groups; i++) {
			num_data_points += raw_data[i].size();
		}
		p_value = dpnht::internal_kruskal_wallis::ComputePValue(num_groups, num_data_points, private_h_statistic, epsilon, mt, p_value_trials);
	}

	// return statistic and p value
	return std::pair<const double, const double>(private_h_statistic, p_value);
}

#endif