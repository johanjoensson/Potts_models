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
		 : index_m(idx), coord_m(), pos_m(pos)
		{
			size_t index = index_m;
			for(size_t i = 0; i < dim; i++)
			{
				size_t offset = 1;
				for(size_t j = 0; j < dim - i - 1; j++ ){
					offset *= size[j];
				}
				coord_m[i] = index / offset;
				index -= coord_m[i]*offset;
			}
		}

		size_t index() const {return index_m;}
		std::array<size_t, dim> coord() const {return coord_m;}
		GSL::Vector pos() const {return pos_m;}
		void set_pos(const GSL::Vector& pos){pos_m = pos;}
};

#endif // SITE_H
