#ifndef SITE_H
#define SITE_H
#include "GSLpp/vector.h"

template<size_t dim>
class Site_t;

namespace std {
	template<size_t dim>
	struct hash<Site_t<dim>>{
		hash() = default;
		size_t operator()(Site_t<dim> &s) const
		{
			return GSL::Vector_hasher_t<double, gsl_vector, std::allocator<double>>()(s.pos()) ^
				std::hash<size_t>()(s.index());
		}
	};
}

template<size_t dim>
class Site_t{
	private:
		size_t index_m;
		std::array<size_t, dim> coord_m;
		GSL::Vector pos_m;
	public:
		Site_t() : index_m(), coord_m(), pos_m() {}
		Site_t(const size_t idx, const GSL::Vector& pos, const std::array<size_t, dim>& size)
		 : index_m(idx), coord_m(calc_coord(idx, size)), pos_m(pos)
		{}

		Site_t(const std::array<size_t, dim>& coords, const GSL::Vector& pos, const std::array<size_t, dim>& size)
		 : index_m(calc_index(coords, size)), coord_m(coords), pos_m(pos)
		{}

		std::array<size_t, dim> calc_coord(const size_t index, const std::array<size_t, dim>& size)
		{
			size_t tmp = index;
			std::array<size_t, dim> res;
			for(size_t i = dim; i > 0; i--)
			{
				size_t offset = 1;
				for(size_t j = 0; j < i - 1; j++ ){
					offset *= size[j];
				}
				res[i - 1] = tmp / offset;
				tmp -= res[i - 1]*offset;
			}
			return res;
		}

		size_t calc_index(const std::array<size_t, dim>& coords, const std::array<size_t, dim>& size)
		{
			size_t res = 0;
			for(size_t i = dim; i > 0; i--)
			{
				size_t offset = 1;
				for(size_t j = 0; j < i - 1; j++ ){
					offset *= size[j];
				}
				res += offset * coords[i - 1];
			}
			return res;
		}

		size_t index() const {return index_m;}
		std::array<size_t, dim> coord() const {return coord_m;}
		GSL::Vector pos() const {return pos_m;}
		void set_pos(const GSL::Vector& pos){pos_m = pos;}

		bool operator==(const Site_t<dim>& s) const
		{
			return (index() == s.index() && coord() == s.coord() && pos() == s.pos());
		}

		bool operator!=(const Site_t<dim>& s) const
		{
			return !(*this == s);
		}
};

template<int dim>
struct Site_t_hasher
{
	size_t operator()(const Site_t<dim>& site) const
	{
		size_t array_hash = 0;
		for(auto val : site.coord()){
			array_hash ^ val;
		}
		return ((std::hash<size_t>()(site.index())^
		(array_hash << 1)) >> 1)^
		(GSL::Vector_hasher()(site.pos()) << 1);
	}
};
#endif // SITE_H
