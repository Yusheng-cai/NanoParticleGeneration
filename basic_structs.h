#pragma once
#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Eigen/Eigen"
struct Atom
{
	float mass = 0.0, charge = 0.0;
	Eigen::Vector3d pos; //kludge to do rotations more easily
	float vx, vy, vz;
	std::string name, type, resname;
	int index, resnr, chargegrp;
	int neighbor_counter = 0;
	Atom();
	Atom(std::string iname, float xi, float yi, float zi);
};

struct Bond
{
	// 1    0.15290   224262.4
	float b0 = 0.15290;
	float kb = 224262.4;
	float a = 0;
	int function = 1;
	int i, j;
	bool has_params = 1;
};
struct Pair
{
	int i;
	int j;
	int function;
};
struct Angle
{
	int function = 1;
	int i, j, k;
	float angle;
	float kb;
};

struct Dihedral
{
	int function = 9;
	int i, j, k, l;
	float angle;
	float kb;
	float multiplicity;
};

struct Ligand
{
	std::vector<Atom> atoms;
	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;
	std::vector<Pair> exclusions;
	std::string name;
	int nrexcl;
	int index;
	Eigen::Vector3d axis;
	void add(Atom a1);
	void shift(float x, float y, float z);
	void set_axis();
	void rotate_to(Eigen::Vector3d target_axis);
	void add_bond(int i, int j);
	void add_angle(int i, int j, int k);
	void add_dihedral(int i, int j, int k, int l);
	void generate_exclusions(int nr_excl);
	void get_exclusions_from(int idx, int nr_excl, std::vector< std::vector<int> >& pair_list);
	void get_neighboring_idx(int idx, std::vector<int>& idx_list, std::vector<int>& dist_list);
};

struct Residue
{
	std::vector<Atom> atoms;
	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;
	std::vector<Pair> exclusions;
	std::string name;
	int index;
	int offset;
	int nrexcl;
	
};


struct Box
{
	float bx, by, bz;
	std::vector<Atom> atoms;
	std::vector<Residue> res;
	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;
	std::vector<Pair> exclusions;
};

void rotate_to(std::vector<Atom>& atoms, Eigen::Vector3d axis, Eigen::Vector3d target_axis);
Ligand load_ligand(std::string topname, std::string groname);
Residue lig2res(Ligand l);
int read_gro_np(std::string filename, std::vector<Atom>& atoms);
int read_gro_np2(std::string filename, std::vector<Atom>& atoms, std::vector<Atom>& sulfurs);
Ligand make_c12_ligand();