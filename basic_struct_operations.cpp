#include "basic_structs.h";
Atom::Atom(std::string iname, float xi, float yi, float zi)
{
	name = iname;
	pos[0] = xi;
	pos[1] = yi;
	pos[2] = zi;
	vx = 0.0f;
	vy = 0.0f;
	vz = 0.0f;
	index = 0;
}

Atom::Atom()
{
	name = "";
	pos[0] = 0.0f; pos[1] = 0.0f; pos[2] = 0.0f;
	vx = 0.0f; vy = 0.0f; vz = 0.0f;
	index = 0;
}

void Ligand::add(Atom a1)
{
	int size = atoms.size();
	a1.index=size;
	atoms.push_back(a1);
	return;
}
void Ligand::shift(float x, float y, float z)
{
	for(int i = 0; i < atoms.size(); i++)
	{
		atoms[i].pos[0] += x;
		atoms[i].pos[1] += y;
		atoms[i].pos[2] += z;
	}
	return;
}
void Ligand::set_axis()
{
	int i2 = atoms.size()-1;
	float x1 = atoms[0].pos[0]; float y1 = atoms[0].pos[1]; float z1 = atoms[0].pos[2];
	float x2 = atoms[i2].pos[0]; float y2 = atoms[i2].pos[1]; float z2 = atoms[i2].pos[2];
	axis[0] = x2-x1;
	axis[1] = y2-y1;
	axis[2] = z2-z1;
	axis.normalize();
	return;
}

Eigen::Matrix3d rotateAlign( Eigen::Vector3d u1, Eigen::Vector3d u2)
{
    Eigen::Vector3d axis = u1.cross(u2);
//	std::cout << u1 << std::endl << u2 << std::endl;
    float cosA = u1.dot(u2);
	float k = 1.0f / (1.0f + cosA);

    Eigen::Matrix3d result;
	result <<
	(axis[0] * axis[0] * k) + cosA,    (axis[1] * axis[0] * k) - axis[2], (axis[2] * axis[0] * k) + axis[1],
	(axis[0] * axis[1] * k) + axis[2], (axis[1] * axis[1] * k) + cosA,    (axis[2] * axis[1] * k) - axis[0],
	(axis[0] * axis[2] * k) - axis[1], (axis[1] * axis[2] * k) + axis[0], (axis[2] * axis[2] * k) + cosA;

    return result;
}
void Ligand::rotate_to(Eigen::Vector3d target_axis)
{
	target_axis.normalize();
	axis.normalize();
	Eigen::Matrix3d R;
	R = rotateAlign(axis, target_axis);
	//std::cout << R;
	for(int i = 0; i < atoms.size(); i++)
	{
		atoms[i].pos = R*atoms[i].pos;
	}
	
	return;
}

void rotate_to(std::vector<Atom>& atoms, Eigen::Vector3d axis, Eigen::Vector3d target_axis)
{
	target_axis.normalize();
	axis.normalize();
	Eigen::Matrix3d R;
	R = rotateAlign(axis, target_axis);
	//std::cout << R;
	for(int i = 0; i < atoms.size(); i++)
	{
		atoms[i].pos = R*atoms[i].pos;
	}	
	return;
}

void Ligand::add_bond(int i, int j)
{
	Bond b1;
	b1.i = i;
	b1.j = j;
	b1.function = 1;
	bonds.push_back(b1);
}
void Ligand::add_angle(int i, int j, int k)
{
	Angle a1;
	a1.i = i;
	a1.j = j;
	a1.k = k;
	a1.function = 1;
	angles.push_back(a1);	
}
void Ligand::add_dihedral(int i, int j, int k, int l)
{
	Dihedral d1;
	d1.i = i;
	d1.j = j;
	d1.k = k;
	d1.l = l;
	d1.function = 9;
	dihedrals.push_back(d1);		
}
void Ligand::get_neighboring_idx(int it, std::vector<int>& idx_list, std::vector<int>& dist_list)
{
	int idx = idx_list[it];
	int idx1, idx2;
	for(int i = 0; i < bonds.size(); i++)
	{
		idx1 = bonds[i].i;
		idx2 = bonds[i].j;
		if(idx == idx1)
		{
			idx_list.push_back(idx2);
			dist_list.push_back(dist_list[it] + 1);
		}
		else if(idx == idx2)
		{
			idx_list.push_back(idx1);
			dist_list.push_back(dist_list[it] + 1);
		}
	}
	return;
}

void Ligand::get_exclusions_from(int idx, int nr_excl, std::vector< std::vector<int> >& pair_list)
{
	std::vector<int> idx_list;
	std::vector<int> dist_list;
	int dist_it = 0;
	idx_list.push_back(idx);
	dist_list.push_back(dist_it);
	for(int i = 0; i < idx_list.size(); i++)
	{
		get_neighboring_idx(i,idx_list,dist_list);
		if(dist_list[i] >= nr_excl) break;
	}
	for(int i = 0; i < idx_list.size(); i++) //remove duplicates with higher distance values
	{
		for(int j = i+1; j < idx_list.size(); j++)
		{
			if(idx_list[i] == idx_list[j])
			{
				idx_list.erase(idx_list.begin() + j);
				dist_list.erase(dist_list.begin() + j);
				j--;
			}
		}
	}
	for(int i = 0; i < idx_list.size(); i++)
	{
		std::vector<int> pair(2,0);
		pair[0] = idx;
		pair[1] = idx_list[i];
		pair_list.push_back(pair);
	}
	return;
}

void Ligand::generate_exclusions(int nr_excl)
{
	std::vector< std::vector<int> > pair_list;
	std::vector<int> distance(atoms.size(), 0);
	for(int i = 0; i < atoms.size(); i++)
	{
		get_exclusions_from(i, nr_excl, pair_list);
	}

	for(int i = 0; i < pair_list.size(); i++)
	for(int j = i+1; j < pair_list.size(); j++)
	{
		int i11 = pair_list[i][0], i12 = pair_list[i][1], i21 = pair_list[j][0], i22 = pair_list[j][1];
		if( (i11 == i21 && i12 == i22) || (i11 == i22 && i12 == i21) )
		{
			pair_list.erase(pair_list.begin() + j);
			j--;
		}
	}
	for(int i = 0; i < pair_list.size(); i++)
	{
		Pair p1;
		p1.i = pair_list[i][0];
		p1.j = pair_list[i][1];
		exclusions.push_back(p1);
	}
	return;
}

Residue lig2res(Ligand l)
{
	Residue r;
	r.atoms = l.atoms;
	r.bonds = l.bonds;
	r.angles = l.angles;
	r.dihedrals = l.dihedrals;
	r.exclusions = l.exclusions;
	r.index = l.index;
	r.name = l.name;
	r.nrexcl = l.nrexcl;
	return r;
}