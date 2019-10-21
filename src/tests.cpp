#include "tests.h"
#include <cstdio>
#include <time.h>
#include <iostream>
#include "safio.h"
#include "potentials.h"

void test_cache()
{
    std::cout << "Testing potential cache speeds" << std::endl;
	double v_tot = 0;
	double v_r_tot = 0;
	double r;
	double r_min = settings.DR_MIN_TAB*8.0 / 1e8;
	double dt_cache;
	double dt_compu;
	int n = 0;
	//Testing stuff goes here.

	clock_t timer = clock();
	//Testing Vr_r calls.
	for (int i = 1; i <= 1e8; i++)
	{
		r = 1 + i * r_min;
		v_tot += Vr_r(r, 1);
		v_r_tot += dVr_dr(r, 1);
		n++;
	}
	dt_cache = ((double)clock() - timer) / CLOCKS_PER_SEC;
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
	dt_compu = ((double)clock() - timer) / CLOCKS_PER_SEC;
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
    for (auto x : lattice.cell_map)
    {
        Cell cell = x.second;
        //Find first un-filled cell.
        if (cell.num == 0) continue;

        //Could use the second, but lets just pull from map to be sure.
        Site &a = lattice.cell_map[x.first].sites[0];
        Site &b = new_lat.cell_map[x.first].sites[0];
        std::cout << &a << " " << &b << std::endl;

        //r_0 is not copied, so we need to test r.
        a.r[0] = 5;
        b.r[0] = 0;
        //If the copy is actually a full copy, this will not be the same.
        counter += a.r[0] - b.r[0];
        break;
    }
    std::cout << "Did it copy?: " << (counter!=0) << std::endl;

    std::cout << "Testing lattice copying speeds" << std::endl;
    clock_t timer = clock();
    int n = 1e4;
    //This seems to work, but leaks memory like crazy...
    for (int i = 0; i < n; i++)
    {
        new_lat = lattice;
    }
    double dt = ((double)clock() - timer) / CLOCKS_PER_SEC;
    dt /= n;
    //Convert to ms;
    dt *= 1000;
    std::cout << " execution time: " << dt << "ms"<< std::endl;
}