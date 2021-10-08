#include "place_sulfurs.h"
int get_grid_size(float& xmin, float& ymin, float& zmin, float& xmax, float& ymax, float& zmax, std::vector<Atom>& atoms)
{
	if(atoms.size() == 0) return -1;
	xmin=atoms[0].pos[0]; ymin=atoms[0].pos[1]; zmin=atoms[0].pos[2];
	xmax=atoms[0].pos[0]; ymax=atoms[0].pos[1]; zmax=atoms[0].pos[2];
	for(int i = 0; i < atoms.size(); i++)
	{
		if(atoms[i].pos[0] < xmin) xmin = atoms[i].pos[0];
		if(atoms[i].pos[0] > xmax) xmax = atoms[i].pos[0];
		
		if(atoms[i].pos[1] < ymin) ymin = atoms[i].pos[1];
		if(atoms[i].pos[1] > ymax) ymax = atoms[i].pos[1];
		
		if(atoms[i].pos[2] < zmin) zmin = atoms[i].pos[2];
		if(atoms[i].pos[2] > zmax) zmax = atoms[i].pos[2];
	}
	xmin -= 0.5;
	xmax += 0.5;
	ymin -= 0.5;
	ymax += 0.5;	
	zmin -= 0.5;
	zmax += 0.5;
	
	return 0;
}

float calc_potential(float x, float y, float z, std::vector<Atom> atoms)
{
	float D = 36.664f,  a = 14.7f, ro = 0.265f;
	float eval = 0.0f, r = 0.0f;
	for(int i = 0; i < atoms.size(); i++)
	{
		if(atoms[i].name == "AU")
		{
			r = (atoms[i].pos[0] - x)*(atoms[i].pos[0] - x) + (atoms[i].pos[1] - y)*(atoms[i].pos[1] - y) + (atoms[i].pos[2] - z)*(atoms[i].pos[2] - z);
			if(r < 1.0f)
			{
				r = sqrt(r);
				eval += D*( (1.0f-std::exp( -a*(r-ro) )) * (1.0f-std::exp( -a*(r-ro) )) ) - D;
			}
		}
	}
	return eval;
}

float calc_potential2(float x, float y, float z, Atom atom)
{
	float D = 36.664f,  a = 14.7f, ro = 0.265f, r = 0.0f;
	r = (atom.pos[0] - x)*(atom.pos[0] - x) + (atom.pos[1] - y)*(atom.pos[1] - y) + (atom.pos[2] - z)*(atom.pos[2] - z);
	r = sqrt(r);
	return D*( (1.0f-std::exp( -a*(r-ro) )) * (1.0f-std::exp( -a*(r-ro) )) ) - D;
}

int calc_potential_v2(std::vector<std::vector<std::vector<float> > >& p_grid, std::vector<Atom>& atoms, float gsx, float gsy, float gsz, float xmin, float ymin, float zmin )
{
	int nx = p_grid.size(); if(nx == 0) return -1;
	int ny = p_grid[0].size(); if(ny == 0) return -1;
	int nz = p_grid[0][0].size(); if(nz == 0) return -1;
	for(int r = 0; r < atoms.size(); r++)
	{
			float xo = atoms[r].pos[0];
			float yo = atoms[r].pos[1];
			float zo = atoms[r].pos[2];
			
			int imin = std::max( (int)floor((xo-xmin)/gsx) - (int)ceil(0.5f/gsx), 0 );
			int jmin = std::max( (int)floor((yo-ymin)/gsy) - (int)ceil(0.5f/gsy), 0 );
			int kmin = std::max( (int)floor((zo-zmin)/gsz) - (int)ceil(0.5f/gsz), 0 );
			
			int imax = std::min( (int)ceil((xo-xmin)/gsx) + (int)ceil(0.5f/gsx), nx-1 );
			int jmax = std::min( (int)ceil((yo-ymin)/gsy) + (int)ceil(0.5f/gsy), ny-1 );
			int kmax = std::min( (int)ceil((zo-zmin)/gsz) + (int)ceil(0.5f/gsz), nz-1 );
			
			float x_temp, y_temp, z_temp;
			for(int i = imin; i <= imax; i++)
			{
				x_temp = gsx*(float)i + xmin;
				for(int j = jmin; j <= jmax; j++)
				{
					y_temp = gsy*(float)j + ymin;
					for(int k = kmin; k <= kmax; k++)
					{
						z_temp = gsz*(float)k + zmin;
						p_grid[i][j][k] += calc_potential2(x_temp, y_temp, z_temp, atoms[r]);
					}
				}
			}
	}
	
	
	return 0;
}


void dump_array(std::vector<std::vector<std::vector<float> > >& p_grid, float gsx, float gsy, float gsz)
{
	std::ofstream ofile("grid_dump.txt");
	for(int i = 0; i < p_grid.size(); i++)
		for(int j = 0; j < p_grid[0].size(); j++)
			for(int k = 0; k < p_grid[0][0].size(); k++)
			{
				ofile << i*gsx << "     " << j*gsy << "     " << k*gsz << "     " << p_grid[i][j][k] << "\n";
			}
	ofile.close();
	return;
}

float get_min_energy_site(std::vector<std::vector<std::vector<float> > >& p_grid, int& ix, int& iy, int& iz)
{
	float e_min = p_grid[0][0][0];
	for(int i = 0; i < p_grid.size(); i++)
		for(int j = 0; j < p_grid[0].size(); j++)
			for(int k = 0; k < p_grid[0][0].size(); k++)
			{
				if(p_grid[i][j][k] < e_min)
				{
					e_min = p_grid[i][j][k];
					ix = i;
					iy = j;
					iz = k;
				}
			}
	return e_min;
}

static inline float get_dist(float a, float b, float c, float x, float y, float z)
{
	return sqrt((a-x)*(a-x) + (b-y)*(b-y) + (c-z)*(c-z));
}

bool get_low_energy_site_near(std::vector<std::vector<std::vector<float> > >& p_grid, int& ix, int& iy, int& iz, int xo, int yo, int zo, float gsx, float gsy, float gsz, float dthresh, float ethresh)
{
	int imin = std::max( xo - (int)ceil(dthresh/gsx), 0 );
	int jmin = std::max( yo - (int)ceil(dthresh/gsy), 0 );
	int kmin = std::max( zo - (int)ceil(dthresh/gsz), 0 );
	
	int imax = std::min( xo + (int)ceil(dthresh/gsx), (int)p_grid.size()-1 );
	int jmax = std::min( yo + (int)ceil(dthresh/gsy), (int)p_grid[0].size()-1 );
	int kmax = std::min( zo + (int)ceil(dthresh/gsz), (int)p_grid[0][0].size()-1 );
	bool hit = 0;
	float e_min = 1e10;
	
	for(int i = imin; i <= imax; i++)
		for(int j = jmin; j <= jmax; j++)
			for(int k = kmin; k <= kmax; k++)
			{
				if(p_grid[i][j][k] < ethresh && get_dist(i*gsx,j*gsy,k*gsz,xo*gsx,yo*gsy,zo*gsz) < dthresh)
				{
					if(p_grid[i][j][k] < e_min)
					{
						e_min = p_grid[i][j][k];
						ix = i;
						iy = j;
						iz = k;
						hit = 1;
					}
				}
			}
	return hit;
}

void dump_array(std::string oname, std::vector<std::vector<std::vector<float> > >& p_grid, float xmin, float ymin, float zmin, float gsx, float gsy, float gsz)
{
	std::ofstream ofile(oname);
	for(int i = 0; i < p_grid.size(); i++)
		for(int j = 0; j < p_grid[0].size(); j++)
			for(int k = 0; k < p_grid[0][0].size(); k++)
			{
				ofile << i*gsx + xmin << "     " << j*gsy + ymin << "     " << k*gsz + zmin << "     " << p_grid[i][j][k] << "\n";
			}
	ofile.close();
	return;
}


float ss_lj2(float r2)
{
	float sigma2 = 0.425*0.425;
	//float eps = 0.39743*4.184;
	if(r2 < sigma2) return 500.0f;
	return 0.0f;
	//return 4.0f*eps*(std::pow(sigma2/r2, 6) - std::pow(sigma2/r2, 3));
}

void add_sulfur(std::vector<std::vector<std::vector<float> > >& p_grid, std::vector<Atom>& sulfurs, float xmin, float ymin, float zmin, float gsx, float gsy, float gsz, int ix, int iy, int iz)
{
	float x = xmin + ix*gsx;
	float y = ymin + iy*gsy;
	float z = zmin + iz*gsz;
	
	Atom a1;
	a1.pos[0] = x; a1.pos[1] = y; a1.pos[2] = z;
	a1.name = "S";
	sulfurs.push_back(a1);
	
	int x_span = ceil(1.0f/gsx);
	int y_span = ceil(1.0f/gsy);
	int z_span = ceil(1.0f/gsz);
	
	int ix_min = std::max(ix-x_span, 0);
	int ix_max = std::min(ix+x_span, (int)p_grid.size() - 1);
	int iy_min = std::max(iy-y_span, 0);
	int iy_max = std::min(iy+y_span, (int)p_grid[0].size() - 1);
	int iz_min = std::max(iz-z_span, 0);
	int iz_max = std::min(iz+z_span, (int)p_grid[0][0].size() - 1);
	
	float x_temp, y_temp, z_temp, r2;
	for(int i = ix_min; i <= ix_max; i++)
	{
		x_temp = i*gsx + xmin;
		for(int j = iy_min; j <= iy_max; j++)
		{
			y_temp = j*gsy + ymin;
			for(int k = iz_min; k <= iz_max; k++)
			{
				z_temp = k*gsz + zmin;
				r2 = (x-x_temp)*(x-x_temp) + (y-y_temp)*(y-y_temp) + (z-z_temp)*(z-z_temp);
				p_grid[i][j][k] += ss_lj2(r2);
			}
		}
	}
	
	
	return;
}

void add_sulfurs(std::vector<std::vector<std::vector<float> > >& p_grid, std::vector<Atom>& sulfurs, float xmin, float ymin, float zmin, float gsx, float gsy, float gsz, float e_cut)
{
	std::string name = "grid2.txt";
	float x, y, z;
	int ix, iy, iz;
	int xo, yo, zo;
	float e_min;
	int sulfur_iterator = 0;
	float energy_thresh = -e_cut;
	float dist_thresh = 0.7f;
	while(true)
	{
		e_min = get_min_energy_site(p_grid, ix, iy, iz);
		if(e_min > energy_thresh) break;
		add_sulfur(p_grid, sulfurs, xmin, ymin, zmin, gsx, gsy, gsz, ix, iy, iz);
		xo = ix; yo = iy; zo = iz;
		while(true)
		{
			if(get_low_energy_site_near(p_grid, ix, iy, iz, xo, yo, zo, gsx, gsy, gsz, dist_thresh, energy_thresh))
			{
				add_sulfur(p_grid, sulfurs, xmin, ymin, zmin, gsx, gsy, gsz, ix, iy, iz);
			}
			else
			{
				if(++sulfur_iterator >= sulfurs.size()) break;
				xo = round((sulfurs[sulfur_iterator].pos[0] - xmin)/gsx);
				yo = round((sulfurs[sulfur_iterator].pos[1] - ymin)/gsy);
				zo = round((sulfurs[sulfur_iterator].pos[2] - zmin)/gsz);
			}
		}
	}
	//dump_array(name, p_grid, xmin, ymin, zmin, gsx, gsy, gsz);
	std::cout << "Number of sulfurs placed = " << sulfurs.size() << std::endl;
	return;
}


int generate_potential_grid(std::vector<Atom>& atoms, std::vector<Atom>& sulfurs, float e_cut)
{
	int nx, ny, nz;
	float xmin, ymin, zmin, xmax, ymax, zmax, gsx, gsy, gsz;


	if(get_grid_size(xmin, ymin, zmin, xmax, ymax, zmax, atoms) == -1) return -1;
	nx = round(-xmin+xmax / 0.02);
	ny = round(-ymin+ymax / 0.02);
	nz = round(-zmin+zmax / 0.02);
	std::cout << "Number of gridpoints: " << nx  << "   " << ny << "   " << nz << std::endl;
	std::vector<std::vector<std::vector<float> > > potential_grid(nx, std::vector< std::vector<float> >(ny, std::vector<float>(nz, 0.0f)) );
	gsx = (xmax - xmin)/((float)(nx - 1));
	gsy = (ymax - ymin)/((float)(ny - 1));
	gsz = (zmax - zmin)/((float)(nz - 1));
	calc_potential_v2(potential_grid, atoms, gsx, gsy, gsz, xmin, ymin, zmin);
	add_sulfurs(potential_grid, sulfurs, xmin, ymin, zmin, gsx, gsy, gsz, e_cut);
	return 0;
}