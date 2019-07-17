#include <bitset>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include "lattice.h"
#include "crystal.h"
#include "potts.h"
#include "GSLpp/error.h"

struct Bitmap_header{
	uint16_t bm;
	uint32_t file_size;
	std::array<uint16_t, 2> reserved;
	uint32_t offset;
	Bitmap_header()
		: bm(0x4D42), file_size(0), reserved{0}, offset(0)
	{}

	uint32_t byte_size()
	{
		return 14;
	}

	void write_to_file(std::ofstream& file)
	{
		file.write(reinterpret_cast<char*>(&bm), 2);
		file.write(reinterpret_cast<char*>(&file_size), 4);
		file.write(reinterpret_cast<char*>(&reserved), static_cast<std::streamsize>(2*reserved.size()));
		file.write(reinterpret_cast<char*>(&offset), 4);
	}
};

std::ostream& operator<<(std::ostream& os, const Bitmap_header& bmh)
{
	os << bmh.bm << " " << bmh.file_size << " ";
	for(auto val : bmh.reserved)
		os << val << " ";
        os << bmh.offset;
	return os;
}

struct DIB_header{
	uint32_t header_size;
	int32_t width, height;
	uint16_t planes, bit_per_pixel;
	uint32_t compression, bitmap_size;
	int32_t h_resolution, v_resolution;
	uint32_t colors, important_colors;
	uint32_t red_bitmask, green_bitmask, blue_bitmask,
		alpha_bitmask, color_space;
	std::array<uint32_t, 16> reserved;
	DIB_header()
		: header_size(124), width(0), height(0), planes(1), bit_per_pixel(32),
		  compression(3), bitmap_size(0), h_resolution(2835),
		  v_resolution(2835), colors(0), important_colors(0),
		  red_bitmask(0x000000FF), green_bitmask(0x0000FF00),
		  blue_bitmask(0x00FF0000), alpha_bitmask(0xFF000000),
		  color_space(0x73524742), reserved{0}
	{}

	uint32_t byte_size()
	{
		return 124;
	}

	void write_to_file(std::ofstream& file)
	{
		file.write(reinterpret_cast<char*>(&header_size), 4);
		file.write(reinterpret_cast<char*>(&width), 4);
		file.write(reinterpret_cast<char*>(&height), 4);
		file.write(reinterpret_cast<char*>(&planes), 2);
		file.write(reinterpret_cast<char*>(&bit_per_pixel), 2);
		file.write(reinterpret_cast<char*>(&compression), 4);
		file.write(reinterpret_cast<char*>(&bitmap_size), 4);
		file.write(reinterpret_cast<char*>(&h_resolution), 4);
		file.write(reinterpret_cast<char*>(&v_resolution), 4);
		file.write(reinterpret_cast<char*>(&colors), 4);
		file.write(reinterpret_cast<char*>(&important_colors), 4);
		file.write(reinterpret_cast<char*>(&red_bitmask), 4);
		file.write(reinterpret_cast<char*>(&green_bitmask), 4);
		file.write(reinterpret_cast<char*>(&blue_bitmask), 4);
		file.write(reinterpret_cast<char*>(&alpha_bitmask), 4);
		file.write(reinterpret_cast<char*>(&color_space), 4);
		file.write(reinterpret_cast<char*>(&reserved), 4*static_cast<std::streamsize>(reserved.size()));
	}
};

struct Pixel_data{
	std::array<uint8_t, 4> color;

	// Bitmap files store color in BGR order
	uint8_t& red(){ return color[0];}
	uint8_t& green(){ return color[1];}
	uint8_t& blue(){ return color[2];}
	uint8_t& alpha(){ return color[3];}

	uint32_t byte_size(){ return 4;}

	void write_to_file(std::ofstream& file)
	{
		file.write(reinterpret_cast<char*>(&color), 4);
	}
};

std::ostream& operator<<(std::ostream& os, const DIB_header& dibh)
{
	os << dibh.header_size << " " << dibh.width << " " << dibh.height << " ";
	os << dibh.planes << " " << dibh.bit_per_pixel << " ";
	os << dibh.compression << " " << dibh.bitmap_size << " ";
	os << dibh.h_resolution << " " << dibh.v_resolution << " ";
	os << dibh.colors << " " << dibh.important_colors << " ";
	os << dibh.red_bitmask << " ";
	os << dibh.green_bitmask << " " << dibh.blue_bitmask << " ";
	os << dibh.alpha_bitmask;
	return os;
}

void bitmap_print(const std::vector<Pixel_data>& pixels, const std::string filename, const std::array<size_t, 2>& size)
{
	Bitmap_header bmh;
	DIB_header dibh;
	bmh.offset = bmh.byte_size() + dibh.byte_size();
	dibh.width = static_cast<int32_t>(size[0]);
	dibh.height = static_cast<int32_t>(size[1]);
	dibh.bitmap_size = static_cast<uint32_t>(pixels.size()) * Pixel_data().byte_size();
	bmh.file_size = bmh.offset  + dibh.bitmap_size;

	std::ofstream file(filename, std::ios::out | std::ios::binary);
	bmh.write_to_file(file);
	dibh.write_to_file(file);
	for(Pixel_data pixel : pixels){
		pixel.write_to_file(file);
	}

	file.close();
}

std::vector<Pixel_data> create_bitmap_data(const std::vector<uint8_t>& field, const std::array<size_t, 2>& size, const unsigned int q = 2)
{
	std::vector<Pixel_data> res(field.size());

	Pixel_data white{0xFF, 0xFF, 0xFF, 0xFF};
	Pixel_data black{0x00, 0x00, 0x00, 0xFF};
	uint8_t val;

	// Bitmap file data starts with the bottom left pixel
	for(auto row = size[1]; row > 0; row--){
		for(size_t column = 0; column < size[0]; column++){
			val = field[(row - 1)*size[0] + column]*0xFF/(q-1);
			res[(row - 1)*size[0] + column] = {val, val, val, 0xFF};
		}
	}
	return res;
}

int main()
{
	GSL::Error_handler e_handler;
	e_handler.off();

	size_t width = 100, height = 100;
	Lattice_t<2> lat{{{static_cast<double>(width), 0},
		{0, static_cast<double>(height)}}};
	Potts_t<2, 3> potts(lat, {width, height}, true);
	potts.set_beta(.5);

	// One interaction parameter per nearest neighbour shell to consider
	potts.set_interaction_parameters({2.0});

	std::tuple<double, int, double, double> corr;


	std::cout << "Iterations start\n";
	size_t num_iterations = width*height;
	for(size_t it = 0; it < num_iterations; it++){
		potts.update(true);
		if(it % (num_iterations/10) == 0){
			std::cout << "Iteration " << it << ", out of " << num_iterations <<"\n";
			std::cout << "\tAverage energy = " << potts.average_site_energy() << "\n";
			std::cout << "\tMagnetization = " << potts.magnetization() << "\n";
			std::cout << "spin-spin correlators  " <<  potts.measure_spin_correlators().size();
			std::cout << "\n";

			bitmap_print(create_bitmap_data(potts.field(), {width, height}, 3),
				"Ising-" + std::to_string((10*it)/num_iterations) + ".bmp", {width, height});
		}
	}

	std::cout << "Average energy = " << potts.average_site_energy() << "\n";
	std::cout << "Magnetization = " << potts.magnetization() << "\n";
	bitmap_print(create_bitmap_data(potts.field(), {width, height}, 3), "Final.bmp", {width, height});

	return 0;
}
