#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "ising.h"

namespace py = pybind11;
using namespace pybind11::literals;

template<size_t dim>
void declare_ising(py::module& m, const std::string& dim_str)
{
	using Ising_class = Ising_t<dim>;
	std::string Ising_string = "Ising_" + dim_str + "D";
	py::class_<Ising_class>(m, Ising_string.c_str())
		.def(py::init<>())
		.def(py::init<Lattice_t<dim>, std::array<size_t, dim>, bool>(), "l"_a, "s"_a, "periodic"_a = false)
		.def("set_Jij", &Ising_class::set_interaction_parameters)
		.def("set_H", &Ising_class::set_H)
		.def("set_beta", &Ising_class::set_beta)
		.def("add_spin_correlator", (void (Ising_class::*)(const size_t, const double)) &Ising_class::add_spin_correlator)
		.def("add_spin_correlator", (void (Ising_class::*)(const std::array<size_t, dim>&, const double)) &Ising_class::add_spin_correlator)
		.def("update", &Ising_class::update, "cluster"_a = false)
		.def("measure_spin_correlators", &Ising_class::measure_spin_correlators)
		.def("total_energy", &Ising_class::total_energy)
		.def("average_site_energy", &Ising_class::average_site_energy)
		.def("magnetization", &Ising_class::magnetization)
		.def_property_readonly("field", &Ising_class::field);

}


PYBIND11_MODULE(ising, m){
	declare_ising<2>(m, "2");
//	declare_ising<3>(m, "3");
}
