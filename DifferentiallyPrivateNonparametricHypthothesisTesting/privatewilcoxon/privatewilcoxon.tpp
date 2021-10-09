#ifndef DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEWILCOXON_TPP
#define DIFFERENTIALLYPRIVATENONPARAMETRICHYPOTHESISTESTING_PRIVATEWILCOXON_TPP

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/privatewilcoxon/privatewilcoxon.h"


template<std::uniform_random_bit_generator Engine>
double dpnht::internal_wilcoxon::PrivateWilcoxonStatistic(const std::vector<std::pair<double, double>>& dataset, double epsilon, Engine& Eng) {
	int size = dataset.size();
	laplace_distribution  ld{ 0, 2 * size / epsilon };
	double public_w_statistic = PublicWilcoxon(dataset);
	return public_w_statistic + ld(Eng);
}

template<std::uniform_random_bit_generator Engine>
double dpnht::internal_wilcoxon::ComputePValue(const int dataset_size, const double private_w_statistic, const double epsilon, Engine& Eng, const int p_value_trials) {
	int successes = 0;
	std::normal_distribution<double>  nd{ 0, (dataset_size * (dataset_size + 1) * (2 * dataset_size + 1)) / epsilon };
	laplace_distribution  ld{ 0, 2 * dataset_size / epsilon };
	double magnitude_private_w_statistic = std::abs(private_w_statistic);
	for (int i = 0; i < p_value_trials; i++) {
		double private_w_statistic_trial = nd(Eng) + ld(Eng);
		if (std::abs(private_w_statistic_trial) > magnitude_private_w_statistic) {
			successes += 1;
		}
	}
	return static_cast<double>(successes) / p_value_trials;
}

template<std::uniform_random_bit_generator Engine>
std::pair<const double, const double> dpnht::internal_wilcoxon::PrivateWilcoxonHelper(const std::vector<std::pair<double, double>>& dataset, const double epsilon, Engine& Eng, const int p_value_trials) {
	int size = dataset.size();
	double private_w_statistic = PrivateWilcoxonStatistic(dataset, epsilon, Eng);
	double p_value = ComputePValue(size, private_w_statistic, epsilon, Eng, p_value_trials);
	// return statistic and p value
	return std::pair<const double, const double>(private_w_statistic, p_value);
}

template<std::uniform_random_bit_generator Engine>
[[nodiscard]] std::pair<const double, const double> dpnht::PrivateWilcoxon(const std::vector<std::pair<double, double>>& dataset, const double epsilon, Engine& Eng, const int p_value_trials) {
	return dpnht::internal_wilcoxon::PrivateWilcoxonHelper(dataset, epsilon, Eng, p_value_trials);
}

#endif