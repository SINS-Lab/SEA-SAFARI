#include "tests.h"
#include <cstdio>
#include <time.h>
#include <iostream>
#include <random>
#include "safio.h"
#include "temps.h"
#include "potentials.h"
#include "space_math.h"

void test_cache()
{
	std::cout << "Testing potential cache speeds" << std::endl;
	double v_tot = 0;
	double v_r_tot = 0;
	double r;
	double r_min = settings.DR_MIN_TAB * 8.0 / 1e8;
	double dt_cache;
	double dt_compu;
	int n = 0;
	//Testing stuff goes here.

	double timer = clock();
	//Testing Vr_r calls.
	for (int i = 1; i <= 1e8; i++)
	{
		r = 1 + i * r_min;
		v_tot += Vr_r(r, 1);
		v_r_tot += dVr_dr(r, 1);
		n++;
	}
	dt_cache = (clock() - timer) / CLOCKS_PER_SEC;
	dt_cache /= n;
	//Convert to ms;
	dt_cache *= 1000;
	std::cout << v_tot << " " << v_r_tot << " ms per: " << dt_cache << std::endl;

	v_tot = 0;
	v_r_tot = 0;
	n = 0;
	timer = clock();
	//Testing Vr_r_init calls.
	for (int i = 1; i <= 1e8; i++)
	{
		r = 1 + i * r_min;
		v_tot += Vr_r_init(r, 1);
		v_r_tot += dVr_dr_init(r, 1);
		n++;
	}
	dt_compu = (clock() - timer) / CLOCKS_PER_SEC;
	dt_compu /= n;
	//Convert to ms;
	dt_compu *= 1000;
	std::cout << v_tot << " " << v_r_tot << " ms per: " << dt_compu << std::endl;
	std::cout << "compu/cache: " << (dt_compu / dt_cache) << std::endl;
}

void test_lattice_copy(Lattice &lattice)
{
	std::cout << "Testing whether it actually copies" << std::endl;

	//If the lattice is made to copy correctly, this will
	//Result in a new lattice, with the same initial conditions.
	double counter = 0;
	Lattice new_lat = lattice;
	for (auto x : new_lat.cell_map)
	{
		Cell *cell = x.second;
		//Find first un-filled cell.
		if (cell->num == 0)
			continue;

		//Could use the second, but lets just pull from map to be sure.
		Site *a = (lattice.cell_map[x.first]->sites[0]);
		Site *b = (new_lat.cell_map[x.first]->sites[0]);
		std::cout << lattice.cell_map[x.first] << " " << new_lat.cell_map[x.first] << std::endl;
		std::cout << a << " " << b << std::endl;

		//r_0 is not copied, so we need to test r.
		a->r[0] = 5;
		b->r[0] = 0;
		//If the copy is actually a full copy, this will not be the same.
		counter += a->r[0] - b->r[0];
		break;
	}
	std::cout << "Did it copy?: " << (counter != 0) << std::endl;

	std::cout << "Testing lattice copying speeds" << std::endl;
	double timer = clock();
	int n = 1e4;
	int x = 0;
	//This seems to work, but leaks memory like crazy...
	for (int i = 0; i < n; i++)
	{
		Lattice *new_lattice = new Lattice(lattice);
		x += new_lattice->id;
		delete new_lattice;
	}
	double dt = (clock() - timer) / CLOCKS_PER_SEC;
	dt /= n;
	//Convert to ms;
	dt *= 1000;
	std::cout << " execution time: " << dt << "ms per, for " << x << std::endl;

	timer = clock();
	double end = timer + CLOCKS_PER_SEC * 10;
	x = 0;
	while (clock() < end)
	{
		x++;
	}
	std::cout << "waited: " << x << std::endl;
}

void test_lattice_springs(Lattice &lattice)
{
	double timer = clock();
	lattice.init_springs(1);
	double dt = (clock() - timer) / CLOCKS_PER_SEC;
	//Convert to ms;
	dt *= 1000;
	std::cout << "Execution time: " << dt << "ms" << std::endl;

	int min = 1e3;
	int max = -1e3;
	for (auto x : lattice.cell_map)
	{
		Cell *cell = x.second;
		if (cell->num <= 0)
			continue;
		min = std::min(cell->sites[0]->near, min);
		max = std::max(cell->sites[0]->near, max);
	}
	std::cout << "Min: " << min << ", Max: " << max << std::endl;
}

void test_mask(Lattice &lattice)
{
	for (int i = 0; i < lattice.mask.num; i++)
	{
		Point p = lattice.mask.points[i];
		debug_file << p.x << "\t" << p.y << std::endl;
	}
	debug_file << std::endl
			   << std::endl;

	for (double x = settings.XSTART; x <= settings.XSTOP; x += settings.XSTEP)
		for (double y = settings.YSTART; y <= settings.YSTOP; y += settings.YSTEP)
		{
			if (lattice.mask.inside(x, y))
			{
				debug_file << x << "\t" << y << std::endl;
			}
		}
	Lattice lettuce = lattice;
	debug_file << "============================="<< std::endl
			   << std::endl;
	for (double x = settings.XSTART; x <= settings.XSTOP; x += settings.XSTEP)
		for (double y = settings.YSTART; y <= settings.YSTOP; y += settings.YSTEP)
		{
			if (lettuce.mask.inside(x, y))
			{
				debug_file << x << "\t" << y << std::endl;
			}
		}
}

void test_rngs()
{
	std::default_random_engine rng;
	rng.seed(make_seed(2.4563));
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << " " << std::endl;
	rng.seed(make_seed(1.4563));
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << " " << std::endl;
	rng.seed(make_seed(2.4563));
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << " " << std::endl;
	rng.seed(make_seed(1.4563));
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
	std::cout << rng() << std::endl;
}