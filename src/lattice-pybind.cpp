#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "GSLpp/matrix.h"
#include "GSLpp/vector.h"
#include "GSLpp/linalg.h"
#include "site.h"
#include "lattice.h"

namespace py = pybind11;

template<size_t dim>
void declare_site(py::module& m, const std::string& dim_str)
{
	using Site_class = Site_t<dim>;
	std::string Site_string = "Site_" + dim_str + "D";
	py::class_<Site_class>(m, Site_string.c_str())
		.def(py::init<>())
		.def_property_readonly("index", &Site_class::index)
		.def_property_readonly("coord", &Site_class::coord)
		.def_property("pos",(GSL::Vector (Site_class::*)()) &Site_class::pos, &Site_class::set_pos);
}

template<size_t dim, class T = double, class M = GSL::Matrix, class V = GSL::Vector>
void declare_lattice(py::module& m, const std::string& dim_str)
{
	using Lattice_class = Lattice_t<dim, T, M, V>;
	std::string Lattice_string = "Lattice_" + dim_str + "D";
	py::class_<Lattice_class>(m, Lattice_string.c_str())
		.def(py::init<>())
		.def(py::init<M>())
		.def(py::init([](py::buffer b){

			py::buffer_info info = b.request();
			if(info.format != py::format_descriptor<double>::format())
				throw std::runtime_error("Incompatible format: expected a double array!");

			if(info.ndim != 2)
				throw std::runtime_error("Incompatible buffer dimension!");
			M mat = M(info.shape[0], info.shape[1]);
			memcpy(mat.data(), info.ptr, sizeof(double)*mat.tda()*mat.size().second);
			return (new Lattice_class(mat));
		}))
		.def_property_readonly("scale", &Lattice_class::scale)
		.def_property_readonly("lat", &Lattice_class::lat)
		.def_property_readonly("recip_lat", &Lattice_class::recip_lat);

}

PYBIND11_MODULE(lattice, m){
	declare_site<2>(m, "2");
	declare_site<3>(m, "3");
	declare_lattice<2>(m, "2");
	declare_lattice<3>(m, "3");


	py::class_<GSL::Matrix>(m, "Matrix", py::buffer_protocol())
		.def(py::init<size_t, size_t>())
		.def(py::init([](py::buffer b){
			py::buffer_info info = b.request();

			if(info.format != py::format_descriptor<double>::format())
				throw std::runtime_error("Incompatible format: expected a double array!");

			if(info.ndim != 2)
				throw std::runtime_error("Incompatible buffer dimension!");

			GSL::Matrix* new_mat = new GSL::Matrix(info.shape[0], info.shape[1]);
			memcpy(new_mat->data(), info.ptr, sizeof(double)*new_mat->tda()*new_mat->size().second);
			return new_mat;
		}))
		.def("__repr__", &GSL::Matrix::to_string)


		.def_buffer([](GSL::Matrix &mat) -> py::buffer_info{
			return py::buffer_info(
					mat.data(),
					sizeof(double),
					py::format_descriptor<double>::format(),
					2,
					{mat.size().first, mat.tda()},
					{sizeof(double)*mat.tda(), sizeof(double)}
					);
		});
};
