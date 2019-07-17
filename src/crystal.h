#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <vector>
#include <algorithm>
#include <unordered_set>
#include <set>
#include "lattice.h"
#include "site.h"
#include "GSLpp/vector.h"
#include "GSLpp/matrix.h"

template<size_t dim>
using Neighbours =  std::vector<Site_t<dim>>;

template<size_t dim>
class Crystal_t {
	Lattice_t<dim> lat_m;
	std::vector<Site_t<dim>> sites_m;
	std::vector<GSL::Vector> R_m, K_m;
	std::array<size_t, dim> size_m;
public:
	Crystal_t():lat_m(), sites_m(), R_m(), K_m(), size_m(){}
	Crystal_t(const Lattice_t<dim>& lat):lat_m(lat), sites_m(), R_m(), K_m(), size_m(){}
	void add_sites(const std::vector<GSL::Vector>&);
	void set_Rn(const double Rmax);
	void set_Kn(const double Kmax);
	void set_size(const std::array<size_t, dim>& size){size_m = size;}
	std::vector<Neighbours<dim>> calc_nearest_neighbours();
	std::vector<Neighbours<dim>> calc_nearest_neighbours(const size_t n_shells);
	std::vector<std::vector<Neighbours<dim>>> determine_nn_shells(const std::vector<Neighbours<dim>>& nn);

	Lattice_t<dim>& lat(){return lat_m;}
	std::vector<Site_t<dim>>& sites(){return sites_m;}
};


template<size_t dim>
void Crystal_t<dim>::add_sites(const std::vector<GSL::Vector>& positions)
{
	for(size_t i = 0; i < positions.size(); i++){
		sites_m.push_back(Site_t<dim>(i, 1./lat_m.scale() * lat_m.lat() *positions[i], size_m));
	}
}

bool comp_norm(const GSL::Vector& a, const GSL::Vector& b)
{
    return a.norm<double>() < b.norm<double>();
}

template<size_t dim>
bool comp_norm_site(const Site_t<dim>& a, const Site_t<dim>& b)
{
	const GSL::Vector va = a.pos(), vb = b.pos();
	return comp_norm(va, vb);
}

struct Vector_comp_norm{
	bool operator()(const GSL::Vector& a, const GSL::Vector& b){
		return comp_norm(a, b);
	}

};

template<size_t dim>
void Crystal_t<dim>::set_Rn(const double Rmax)
{
	std::vector<int> N(dim, 0);
	GSL::Vector n(dim), tmp;
	GSL::Matrix a, b;
	a.copy(lat_m.lat());
	b.copy(lat_m.recip_lat());
	std::set<GSL::Vector, Vector_comp_norm> res;

	// Calculate limits
	for(size_t i = 0; i < dim; i++){
		N[i] = static_cast<int>(std::ceil(b[i].norm<double>()/(2*M_PI)*Rmax));
		n[i] = 0;
	}

	// temporary vector for storing linear combinations of new and old vector
	std::unordered_set<GSL::Vector, GSL::Vector_hasher_t<double, gsl_vector, std::allocator<double>>> r_tmp;
	for(size_t i = 0; i < dim; i++){
		for(n[i] = -N[i]; n[i] <= N[i]; n[i]++){
			// Add ni * ai to the list of vectors
			R_m.push_back(a*n);
			// add ni*ai to all the vectors already found
			for(GSL::Vector v : R_m){
				r_tmp.insert(a*n + v);
			}
		}
		n[i] = 0;
	}
	R_m.assign(r_tmp.begin(), r_tmp.end());
	std::sort(R_m.begin(), R_m.end(), comp_norm);
}

template<size_t dim>
void Crystal_t<dim>::set_Kn(const double Kmax)
{
	std::vector<int> N(dim, 0);
	GSL::Vector n(dim), tmp;
	GSL::Matrix a, b;
	a.copy(lat_m.lat());
	b.copy(lat_m.recip_lat());
	std::set<GSL::Vector, Vector_comp_norm> res;

	// Calculate limits
	for(size_t i = 0; i < dim; i++){
		N[i] = static_cast<int>(std::ceil(a[i].norm<double>()/(2*M_PI)*Kmax));
		n[i] = 0;
	}

	// temporary vector for storing linear combinations of new and old vector
	std::unordered_set<GSL::Vector, GSL::Vector_hasher> k_tmp;
	for(size_t i = 0; i < dim; i++){
		for(n[i] = -N[i]; n[i] <= N[i]; n[i]++){
			// Add ni * ai to the list of vectors
			K_m.push_back(b*n);
			// add ni*ai to all the vectors already found
			for(GSL::Vector v : K_m){
				k_tmp.insert(b*n + v);
			}
		}
		n[i] = 0;
	}
	K_m.assign(k_tmp.begin(), k_tmp.end());
	std::sort(K_m.begin(), K_m.end(), comp_norm);
}

template<size_t dim>
std::vector<Neighbours<dim>> Crystal_t<dim>::calc_nearest_neighbours()
{
	std::vector<Neighbours<dim>> res(sites_m.size());
	GSL::Vector ri, rj, zero_v(dim);
	for(size_t i = 0; i < sites_m.size(); i++){
		ri = sites_m[i].pos();
		// Loop over all sites inside the cell
		for(size_t j = 0; j < i; j++){
			rj = sites_m[j].pos();
			// Add all lattice vectors
			for(auto R : R_m){
				// Insert all atoms in the system
				res[i].push_back(Site_t<dim>(j, rj + R - ri, size_m));
				res[j].push_back(Site_t<dim>(i, ri - R - rj, size_m));
			}
		}
		std::sort(res[i].begin(), res[i].end(), comp_norm_site<dim>);
	}
	return res;
}

template<>
std::vector<Neighbours<3>> Crystal_t<3>::calc_nearest_neighbours(const size_t n_steps)
{
	return this->calc_nearest_neighbours();
}

template<>
std::vector<Neighbours<2>> Crystal_t<2>::calc_nearest_neighbours(const size_t n_steps)
{
	std::vector<Neighbours<2>> res(sites_m.size());
	GSL::Vector ri, rj, R(2), zerov(2);
	GSL::Matrix a = lat_m.lat();
	size_t j = 0;
	double a0, a1;
	for(size_t i = 0; i < sites_m.size(); i++){
		ri = sites_m[i].pos();
		for(size_t step_h = 0; step_h < n_steps + 1; step_h++){
			for(size_t step_v = 0; step_v < n_steps + 1; step_v++){
				if(step_h == 0 && step_v == 0){
					continue;
				}

				if((i % size_m[0]) + step_h < size_m[0] && (i / size_m[0]) + step_v < size_m[1]){
					j = i + step_v*size_m[0] + step_h;
					R .copy(zerov);
				}else if(R_m.size() != 0){
					j = ((i / size_m[0] + step_v) % size_m[1])*size_m[0] +
						(i + step_h) % size_m[0];
					a0 = static_cast<double>((i % size_m[0] + step_h) / size_m[0]);
					a1 = static_cast<double>((i / size_m[0] + step_v) / size_m[1]);
					R  = a0 * a[0] + a1 * a[1];
				}
				rj = sites_m[j].pos();
				res[i].push_back(Site_t<2>(j, rj + R - ri, size_m));

				if(step_h > 0){
					if((i % size_m[0]) >= step_h   && (i / size_m[0] + step_v) < size_m[1]){
						j = i + step_v*size_m[0] - step_h;
						R.copy(zerov);
					}else if(R_m.size() != 0){
						j = ((i / size_m[0] + step_v) % size_m[1])*size_m[0];
						if(i % size_m[0] < step_h){
							j += (i + (size_m[0] - (step_h % size_m[0]))) % size_m[0];
							a0 = static_cast<double>((step_h - i % size_m[0])/ size_m[0] + 1);
							R = a0*(-a[0]);
						}else{
							j +=  i % size_m[0] - step_h;
							R.copy(zerov);
						}
						a1 = static_cast<double>((i / size_m[0] + step_v) / size_m[1]);
						R += a1*a[1];
					}
					rj = sites_m[j].pos();
					res[i].push_back(Site_t<2>(j, rj + R - ri, size_m));
				}

				if(step_v > 0){
					if((i % size_m[0] + step_h) < size_m[0] && (i / size_m[0]) >= step_v){
						j = i - step_v*size_m[0] + step_h;
						R.copy(zerov);
					}else if(R_m.size() != 0){
						if(i / size_m[0] < step_v){
							j = (i/size_m[0] + (size_m[1] - (step_v % size_m[1])) % size_m[1]  )*size_m[0];
							a1 = static_cast<double>((step_v - i/size_m[1])/ size_m[1] + 1);
							R = a1*(-a[1]);
						}else{
							j = (i/size_m[0] - step_v)*size_m[0];
							R.copy(zerov);
						}
						j += (i + step_h) % size_m[0];
						a0 = static_cast<double>((i % size_m[0] + step_h) / size_m[0]);
						R  +=  a0 * a[0];
					}
					rj = sites_m[j].pos();
					res[i].push_back(Site_t<2>(j, rj + R - ri, size_m));
				}

				if(step_v > 0 && step_h > 0){
					if((i % size_m[0]) >= step_h   && (i / size_m[0] ) >= step_v){
						j = i - step_v*size_m[0] - step_h;
						R.copy(zerov);
					}else if(R_m.size() != 0){
						if(i % size_m[0] < step_h){
							j = (i + size_m[0] - (step_h % size_m[0])) % size_m[0];
							a0 = static_cast<double>((step_h - i % size_m[0])/ size_m[0] + 1);
							R = a0*(-a[0]);
						}else{
							j = (i - step_h) % size_m[0];
							R.copy(zerov);
						}
						if(i/size_m[0] < step_v){
							j += (i/size_m[0] + (size_m[1] - (step_v % size_m[1])) % size_m[1]  )*size_m[0];
							a1 = static_cast<double>((step_v -  i / size_m[0]) / size_m[1] + 1);
							R += a1*(-a[1]);
						}else{
							j += (i/size_m[0] - step_v)*size_m[0];
							R += zerov;
						}
					}
					rj = sites_m[j].pos();
					res[i].push_back(Site_t<2>(j, rj + R - ri, size_m));
				}

			}
		}
		std::sort(res[i].begin(), res[i].end(), comp_norm_site<2>);
	}
	return res;
}

template<size_t dim>
std::vector<std::vector<Neighbours<dim>>> Crystal_t<dim>::determine_nn_shells(const std::vector<Neighbours<dim>>& nn)
{
	std::vector<std::vector<Neighbours<dim>>> res(nn.size());
	size_t shell_idx;
	double old_dist;
	Site_t<dim> tmp;
	for(size_t site_idx = 0; site_idx < nn.size(); site_idx++){
		tmp = nn[site_idx][0];
		old_dist = tmp.pos(). template norm<double>();
		res[site_idx] = std::vector<Neighbours<dim>>(1, Neighbours<dim>());
		shell_idx = 0;
		for(size_t neighbour_idx = 0; neighbour_idx < nn[site_idx].size(); neighbour_idx++){
			tmp = nn[site_idx][neighbour_idx];
			if(std::abs( tmp.pos(). template norm<double>() - old_dist ) > 1e-6){
				old_dist = tmp.pos(). template norm<double>();
				shell_idx++;
				res[site_idx].push_back(Neighbours<dim>());
			}
			res[site_idx][shell_idx].push_back( tmp );
		}
	}

	return res;
}
#endif //CRYSTAL_H
