// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

#ifndef DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEMANNWHITNEY_H
#define DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEMANNWHITNEY_H

// TO DO: make a struct containing epsilon_m, epsilon_u, and delta and pass that as function argument
// (shortens argument list and allows you to pass other things you don't want to recompute)
// TO DO: instead of computing random values and then computing rankings in compute p-value,
// use the random values to randomly assign ranks (e.g., a value of <.5 means group 1 has the next
// rank and a value of >.5 means group 2 has the next rank)

// map should hold pairs of integers with the first integer holding the number in the first group and the second holding the number in the second group


#include <map>
#include <utility>
#include <concepts>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/laplacedistribution/laplacedistribution.h"

namespace dpnht {

	template<std::totally_ordered datatype = double, std::uniform_random_bit_generator Engine = std::mt19937>
	[[nodiscard]] std::pair<const double, const double> PrivateMannWhitney(const std::pair<std::vector<datatype>, std::vector<datatype>>& raw_data, const double epsilon_m, const double epsilon_u, const double delta, Engine& Eng, const int p_value_trials);

	template<std::totally_ordered datatype = double>
	[[nodiscard]] std::pair<const double, const double> PrivateMannWhitney(const std::pair<std::vector<datatype>, std::vector<datatype>>& raw_data, const double epsilon_m, const double epsilon_u, const double delta, const int p_value_trials);


	namespace internal_mann_whitney {

		double PublicMannWhitney(const std::pair<std::vector<double>, std::vector<double>>& rankings_data);

		template<std::uniform_random_bit_generator Engine>
		std::pair<double, double> PrivateMannWhitneyStatistic(const std::pair<std::vector<double>, std::vector<double>>& rankings_data, const double epsilon_m, const double epsilon_U, const double delta, Engine& Eng);


		template<std::uniform_random_bit_generator Engine>
		std::pair<std::vector<double>, std::vector<double>> BuildRandomDatabase(const std::pair<int, int>& group_sizes, Engine& Eng);

		template<std::totally_ordered datatype = double>
		std::pair<std::vector<double>, std::vector<double>> ComputeRankings(const std::pair<std::vector<datatype>, std::vector<datatype>>& raw_data);

		template<std::uniform_random_bit_generator Engine>
		double ComputePValue(const std::pair<int, int>& group_sizes, const double private_u_statistic, const double epsilon_m, const double epsilon_u, const double delta, Engine& Eng, const int p_value_trials);

		template<std::totally_ordered datatype = double, std::uniform_random_bit_generator Engine = std::mt19937>
		std::pair<const double, const double> PrivateMannWhitneyHelper(const std::pair<std::vector<datatype>, std::vector<datatype>>& raw_data, const double epsilon_m, const double epsilon_u, const double delta, Engine& Eng, const int p_value_trials);

	}
}

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/privatemannwhitney/privatemannwhitney.tpp"

#endif
