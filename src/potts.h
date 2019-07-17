#ifndef POTTS_H
#define POTTSG_H

#include <vector>
#include <tuple>
#include <algorithm>
#include <random>

#include <iostream>
#include <iomanip>

#include "crystal.h"
#include "site.h"
#include "lattice.h"
#include "GSLpp/matrix.h"

template<size_t dim, size_t q>
class Potts_t{
	using Site = Site_t<dim>;
	private:
		std::array<size_t, dim> size_m;
		Crystal_t<dim> cr_m;
		std::vector<uint8_t> field_m;
		std::vector<std::vector<Neighbours<dim>>> nn_shells_m;
		std::vector<double> J_m;
		double H_m;
		double beta_m;
		std::vector<std::tuple<size_t, size_t, GSL::Vector>> correlators_m;

		size_t calc_length() const
		{
			size_t length = 1;
			for(auto val : size_m){
				length *= val;
			}
			return length;
		}

		void setup_field()
		{
			size_t length = calc_length();
			field_m = std::vector<uint8_t>(length);
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<uint8_t> dist(0,q - 1);
			for(auto& val : field_m){
				val = static_cast<int8_t>(dist(gen));
			}
		}

		void setup_crystal()
		{
			cr_m.set_size(size_m);
		}

		GSL::Vector calc_coords(const size_t index)
		{
			GSL::Vector res(dim);
			size_t offset = 1;
			size_t tmp_index = index;
			for(size_t i = dim; i > 0; i--)
			{
				offset = 1;
				for(size_t j = 0; j < i - 1; j++ ){
					offset *= size_m[j];
				}
				res[i-1] = static_cast<double>(tmp_index / offset);
				tmp_index -= (tmp_index/offset)*offset;
			}
			return res;
		}

		size_t calc_index(const std::array<size_t, dim>& coords)
		{
			size_t res = 0, offset;
			for(size_t i = dim; i > 0; i--)
			{
				offset = 1;
				for(size_t j = 0; j < i - 1; j++ ){
					offset *= size_m[j];
				}
				res += coords[i]*offset;
			}
			return res;

		}

		void setup_crystal_sites()
		{
			size_t length = calc_length();
			std::vector<GSL::Vector> positions(length);
			GSL::Vector coord(dim);
			for(size_t idx = 0; idx < length; idx++){
				coord = calc_coords(idx);
				positions[idx] = 1./cr_m.lat().scale()*cr_m.lat().lat()*coord;
			}
			cr_m.add_sites(positions);
		}

		void setup_crystal_lattice_vectors(bool periodic)
		{
			if(periodic){
				cr_m.set_Rn(1);
			}else{
				cr_m.set_Rn(0);
			}
		}

		void setup_nearest_neighbour_shells(const size_t n_steps = 1)
		{
			std::vector<Neighbours<dim>> nn = cr_m.calc_nearest_neighbours(n_steps);
			nn_shells_m = std::move(cr_m.determine_nn_shells(nn));
		}

		// Delta function for the interactions
		double site_energy(const size_t index) const
		{
			double energy = 0;
			uint8_t other_spins = 0;
			// Loop over all interaction constants provided
			for(size_t i = 0; i < J_m.size(); i++){
				for(const auto neighbour : nn_shells_m[index][i]){
					if(field_m[neighbour.index()] == field_m[index]){
						other_spins++;
					}
				}
				energy -= J_m[i]/2 * other_spins;
				other_spins = 0;
			}
			// External field aligned with 0th spin state
			if(field_m[index] == 0){
				energy -= H_m;
			}

			return energy;
		}

		uint8_t change_spin(const size_t index) const
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<uint8_t> dist(0,q - 1);
			uint8_t res = field_m[index];
			while(res == field_m[index]){
				res = dist(gen);
			}
			return res;
		}

		void flip_single_spin(const size_t index)
		{
			std::uniform_real_distribution<double> dist_d(0., 1.);
			std::random_device rd;
			std::mt19937 gen(rd());
			uint8_t old_spin = field_m[index];
			double e_site = site_energy(index);
			field_m[index] = change_spin(index);
			double e_trial = site_energy(index);

			if( e_trial > e_site && dist_d(gen) > GSL::exp(-beta_m*(e_trial - e_site)).val){
				field_m[index] = old_spin;
			}
		}

		std::vector<size_t> build_cluster(const size_t index)
		{
			std::vector<size_t> to_treat{index};
			std::unordered_set<size_t> treated;
			std::uniform_real_distribution<double> dist_d(0., 1.);
			std::random_device rd;
			std::mt19937 gen(rd());
			double J;
			size_t i;

			while(to_treat.size() > 0){
				i = to_treat.back();
				to_treat.pop_back();
				for(size_t n_shell = 0; n_shell < J_m.size(); n_shell++){
					J = J_m[n_shell];
					for(auto neighbour : nn_shells_m[i][n_shell]){
						if(J > 0){
							if(field_m[neighbour.index()] == field_m[i] && dist_d(gen) > GSL::exp(-beta_m*J).val && treated.find(neighbour.index()) == treated.end()){
								to_treat.push_back(neighbour.index());
							}
						}else if(J < 0){
							if(field_m[neighbour.index()] != field_m[i] && dist_d(gen) > GSL::exp(beta_m*J).val && treated.find(neighbour.index()) == treated.end()){
								to_treat.push_back(neighbour.index());
							}
						}
					}
					treated.insert(i);
				}
			}
			std::vector<size_t> res(treated.size());
			auto res_it = res.begin();
			for(auto tmp_it = treated.begin(); tmp_it != treated.end(); tmp_it++, res_it++){
				*res_it = *tmp_it;
			}
			return res;
		}

		void flip_spin_cluster(const size_t index)
		{
			auto cluster = build_cluster(index);
			uint8_t new_spin = change_spin(index);
			for(size_t i : cluster){
				field_m[i] = new_spin;
			}
		}

	public:
		Potts_t() : size_m(), cr_m(), field_m(), nn_shells_m(), J_m(), H_m(0), beta_m(), correlators_m() {}
		Potts_t(const Lattice_t<dim>& l, const std::array<size_t, dim> & s, bool periodic = false)
			: size_m(s), cr_m(l), field_m(), nn_shells_m(), J_m(), H_m(0), beta_m(), correlators_m()
		{
			setup_field();
			setup_crystal();
			setup_crystal_sites();
			setup_crystal_lattice_vectors(periodic);
			setup_nearest_neighbour_shells(2);
		}

		void set_interaction_parameters(const std::vector<double>& J)
		{
			J_m = J;
			if(J_m.size() > 1){
				std::cout << "Recalculating nearest neighbours to enable inclusion of (at least) " << J_m.size() << " nearest neighbour shells\n";
				setup_nearest_neighbour_shells(J_m.size()/2 + 1);
			}
		}

		void set_H(const double H){H_m = H;}
		void set_beta(const double beta){beta_m = beta;}

		void add_spin_correlator(const size_t index, const double r_max)
		{
			double r = 0;
			size_t j = 0;
			while(r <= r_max && j < nn_shells_m[index].size() && nn_shells_m[index][j].size() > 0){
				r = nn_shells_m[index][j][0].pos(). template norm<double>();
				for(const auto& site : nn_shells_m[index][j]){
					correlators_m.push_back(std::make_tuple(index, site.index(), site.pos()));
				}
				j++;
			}
		}

		void add_spin_correlator(const std::array<size_t, dim>& i, const double r_max)
		{
			add_spin_correlator(calc_index(i), r_max);
		}

		// void add_spin_correlator(const std::array<size_t, dim>& i, const std::array<size_t, dim>& j)
		// {
		// 	size_t i_index = 0;
		// 	size_t j_index = 0;
		//
		// 	add_spin_correlator(i_index, j_index);
		// }

		std::vector<uint8_t>& field(){return field_m;}

		void update(bool cluster = false)
		{
			size_t length = calc_length();
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<size_t> dist_i(0, length - 1);
			size_t index = dist_i(gen);
			if(cluster){
				flip_spin_cluster(index);
			}else{
				flip_single_spin(index);
			}
		}

		int spin_spin(const size_t i, const size_t j) const
		{
			return field_m[i]*field_m[j];
		}

		std::vector<std::tuple<double, double, double>> measure_spin_correlators()
		{
			std::vector<std::tuple<double, double, double>> res(correlators_m.size());
			size_t i, j;
			unsigned int si, ss, sj;
			double r = 0;
			std::tuple<size_t, size_t, GSL::Vector> corr;
			for(auto index = 0; index < correlators_m.size(); index++){
				corr = correlators_m[index];
				i = std::get<0>(corr);
				j = std::get<1>(corr);
				si = field_m[i];
				sj = field_m[j];
				r = std::get<2>(corr).template norm<double>();
				res[index] = std::make_tuple(r, si, sj);
			}
			return res;
		}

		double total_energy() const
		{
			double res = 0;
			for(size_t i = 0; i < field_m.size(); i++){
				res += site_energy(i);
			}
			return res;
		}

		double average_site_energy() const
		{
			size_t n_sites = calc_length();
			return total_energy()/static_cast<double>(n_sites);
		}

		double magnetization() const
		{
			double length = static_cast<double>(calc_length());
			unsigned int tot_spin = 0;
			for(auto spin : field_m){
				tot_spin += spin;
			}
			return static_cast<double>(tot_spin)/length;
		}
};

template<size_t dim> using Ising_t = Potts_t<dim, 2>;
#endif // POTTS_H
