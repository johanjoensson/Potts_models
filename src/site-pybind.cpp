#include <pybind11/pybind11.h>
#include "site.h"

namespace py = pybind11;

template<size_t dim>
void declare_site(py::module& m, const std::string& dim_str)
{
	using Site_t<dim> = Site_class;
	std::string Site_string = "Site_" + dim_str + "D";
	py::class<Site_class>(m, Site_string.c_str())
		.def(py::init<>())
		.def("index", &Site_class::index)
		.def("coord", &Site_class::coord);
}

PYBIND11_MODULE(lattice, m){
	declare_site<2>(m, "2");
	declare_site<3>(m, "3");
};
