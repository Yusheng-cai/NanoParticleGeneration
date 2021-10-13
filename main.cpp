#include "tools/Assert.h"
#include "tools/InputParser.h"

#include "basic_structs.h"
#include "place_sulfurs.h"


float a = 0.4065f; //gold lattice constant in nm
float x_list[14] = {0.0f,0.0f,0.0f,0.0f,a,a,a,a,0.5f*a,0.5f*a,0.0f,0.5f*a,0.5f*a,a};
float y_list[14] = {0.0f, 0.0f, a, a, 0.0f, 0.0f, a, a, 0.5f*a, 0.0f, 0.5f*a, 0.5f*a, a, 0.5f*a};
float z_list[14] = {0.0f, a, 0.0f, a, 0.0f, a, 0.0f, a, 0.0f, 0.5f*a, 0.5f*a, a, 0.5f*a, 0.5f*a};
void gen_cell_fcc(float x, float y, float z, std::vector<Atom>& atoms, int i, int j, int k) 
{
	Atom atom_template;
	atom_template.name = "AU";
	atom_template.type = "AU";
	atom_template.mass = 196.96;
	atom_template.charge = 0.0f;
	for(int r = 0; r < 14; r++)
	{
		if(i != 0 && x_list[r] == 0.0f) continue;
		if(j != 0 && y_list[r] == 0.0f) continue;
		if(k != 0 && z_list[r] == 0.0f) continue;
		atom_template.pos[0] = x_list[r] + x;
		atom_template.pos[1] = y_list[r] + y;
		atom_template.pos[2] = z_list[r] + z;
		atoms.push_back(atom_template);
	}
	return;
}
void get_com(std::vector<Atom>& atoms, float& x, float& y, float& z)
{
	float x_tot = 0.0f;
	float y_tot = 0.0f;
	float z_tot = 0.0f;
	int counter = 0;
	for(int i = 0; i < atoms.size(); i++)
	{
			x_tot += atoms[i].pos[0];
			y_tot += atoms[i].pos[1];
			z_tot += atoms[i].pos[2];
			counter++;
	}
	
	x = x_tot * 1.0f/(float)counter;
	y = y_tot * 1.0f/(float)counter;
	z = z_tot * 1.0f/(float)counter;
	return;
}
void shift_atoms(std::vector<Atom>& atoms, float x, float y, float z)
{
	int counter = 0;
	for(int i = 0; i < atoms.size(); i++)
	{
		atoms[i].pos[0] += x;
		atoms[i].pos[1] += y;
		atoms[i].pos[2] += z;
	}
	return;
}
void gen_lattice_fcc(float x, float y, float z, std::vector<Atom>& atoms) //generates a lattice of gold atoms with origin at 0,0,0 and extending to x, y, z
{
	int nx = x/a, ny = y/a, nz = z/a;
	float origin_x = 0.0f, origin_y = 0.0f, origin_z = 0.0f;
	for(int i = 0; i < nx; i++)
	for(int j = 0; j < ny; j++)
	for(int k = 0; k < nz; k++){
		origin_x = i*a, origin_y = j*a, origin_z = k*a;
		gen_cell_fcc(origin_x, origin_y, origin_z, atoms, i,j,k);
	}
	return;
}
int count_duplicates(std::vector<Atom>& atoms)
{
	int count = 0;
	for(int i = 0; i < atoms.size(); i++)
	{
		for(int j = i+1; j < atoms.size(); j++)
		{
			if(atoms[i].pos[0] == atoms[j].pos[0] && atoms[i].pos[1] == atoms[j].pos[1] && atoms[i].pos[2] == atoms[j].pos[2]) count++;
		}
	}
	return count;
}
void trim_atoms(std::vector<Atom>& atoms, int number)
{
	int counter = 0;
	int n_delete = atoms.size() - number;
	while( counter < n_delete)
	{
		atoms.pop_back();
		counter++;
	}
	return;
}

void cut_sphere(float r, float x_offset, float y_offset, float z_offset, std::vector<Atom>& atoms)
{
	float x,y,z;
	std::vector<Atom> buffer;
	for(int i = 0; i < atoms.size(); i++)
	{
		x = atoms[i].pos[0]-x_offset;
		y = atoms[i].pos[1]-y_offset;
		z = atoms[i].pos[2]-z_offset;
		if(x*x + y*y + z*z > r*r) continue;
		buffer.push_back(atoms[i]);
	}
	atoms = buffer;
	return;
}

void cut_cylinder(float r, float p1[3], float p2[3], std::vector<Atom>& atoms)
{
	Eigen::Vector3d ra(p1[0],p1[1],p1[2]);
	Eigen::Vector3d rb(p2[0],p2[1],p2[2]);
	Eigen::Vector3d e = rb-ra;
	Eigen::Vector3d rp, num;
	float d;
	std::vector<Atom> buffer;
	for(int i = 0; i < atoms.size(); i++)
	{
		rp[0] = atoms[i].pos[0]; rp[1] = atoms[i].pos[1]; rp[2] = atoms[i].pos[2];
		num = e.cross(rp-ra);
		d = num.norm()/e.norm();
		if( d > r ) continue;
		buffer.push_back(atoms[i]);
	}
	atoms = buffer;
	return;
}

void rotate_atoms_about_com(std::vector<Atom>& atoms, float axis[3], float target_axis[3])
{
	Eigen::Vector3d a(axis[0], axis[1], axis[2]);
	Eigen::Vector3d b(target_axis[0], target_axis[1], target_axis[2]);
	float xcom, ycom, zcom;
	get_com(atoms, xcom, ycom, zcom);
	shift_atoms(atoms, -xcom, -ycom, -zcom);
	rotate_to(atoms, a, b);
	shift_atoms(atoms, xcom, ycom, zcom);
	return;
}

void cut_fcc_above(float z, std::vector<Atom>& atoms)
{
	float x,y,z2;
	std::vector<Atom> buffer;
	for(int i = 0; i < atoms.size(); i++)
	{
		x = atoms[i].pos[0]; y = atoms[i].pos[1]; z2 = atoms[i].pos[2];
		if(z2 > z - x - y) continue;
		buffer.push_back(atoms[i]);
	}
	atoms = buffer;
	return;
}

void cut_fcc_below(float z, std::vector<Atom>& atoms)
{
	float x,y,z2;
	std::vector<Atom> buffer;
	for(int i = 0; i < atoms.size(); i++)
	{
		x = atoms[i].pos[0]; y = atoms[i].pos[1]; z2 = atoms[i].pos[2];
		if(z2 < z - x - y) continue;
		buffer.push_back(atoms[i]);
	}
	atoms = buffer;
	return;
}
void generate_bonds(Residue& gnp_res, float multiplier)
{
	float dx, dy, dz, distance;
	
	for(int i = 0; i < gnp_res.atoms.size(); i++)
	{
		for(int j = i+1; j < gnp_res.atoms.size(); j++)
		{
			dx = gnp_res.atoms[i].pos[0] - gnp_res.atoms[j].pos[0];
			dy = gnp_res.atoms[i].pos[1] - gnp_res.atoms[j].pos[1];
			dz = gnp_res.atoms[i].pos[2] - gnp_res.atoms[j].pos[2];
			distance = sqrt(dx*dx + dy*dy + dz*dz);
			if(distance < multiplier*a) //multiplier is 0.75 when dealing with a generated lattice, nanomodeler tends to have AUS, AUL atoms very far from surface with only core atoms being based on lattice spacing
			{
				Bond b1;
				b1.i = i;
				b1.j = j;
				b1.b0 = distance;
				b1.kb = 124410;
				gnp_res.bonds.push_back(b1);
				gnp_res.atoms[i].neighbor_counter++;
				gnp_res.atoms[j].neighbor_counter++;
			}
		}
	}
	int surf_atom_counter = 0;
	for(int i = 0; i < gnp_res.atoms.size(); i++)
	{
		std::cout << gnp_res.atoms[i].neighbor_counter << std::endl;
		if(gnp_res.atoms[i].neighbor_counter < 12){
			//gnp_res.atoms[i].name = "AUS";
			surf_atom_counter++;
		}
	}
	std::cout << "# Surface atoms = " << surf_atom_counter << std::endl;
	
	return;
}
void dump_atoms_gro(std::string file_header, Box box) //will implicitly assume 1 res per atom
{
	std::vector<Residue> residues = box.res;
	std::ofstream ofile(file_header + ".gro");
	ofile << "Generated by Zachariah Vicars, 10/25/20\n";
	char buffer[256];
	sprintf(buffer, "%5d", box.atoms.size());
	ofile << buffer << "\n";
	float max_x = 0,  max_y = 0, max_z = 0;
	int atom_counter = 1;
	for(int j = 0; j < residues.size(); j++)
	for(int i = 0; i < residues[j].atoms.size(); i++)
	{
		if(residues[j].atoms[i].pos[0] > max_x) max_x = residues[j].atoms[i].pos[0];
		if(residues[j].atoms[i].pos[1] > max_y) max_y = residues[j].atoms[i].pos[1];
		if(residues[j].atoms[i].pos[2] > max_z) max_z = residues[j].atoms[i].pos[2];
		sprintf(buffer, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f", j+1, residues[j].name.c_str(), residues[j].atoms[i].name.c_str(), atom_counter, residues[j].atoms[i].pos[0], residues[j].atoms[i].pos[1], residues[j].atoms[i].pos[2], 0.0f, 0.0f, 0.0f);
		std::string line = buffer;
		ofile << buffer << "\n";
		atom_counter++;
	}
	sprintf(buffer, "%10.5f%10.5f%10.5f\n", max_x, max_y, max_z);
	ofile << buffer;
	ofile.close();
	return;
}

void write_itp(std::string file_header, Box box) //will implicitly assume 1 res per atom
{
	std::ofstream ofile(file_header + ".itp");
	ofile << ";Generated by Zachariah Vicars, 10/30/20\n";
	char buffer[256];
	ofile << "[ moleculetype ]\n";
	ofile << ";name       nrexcl\n";
	ofile << "NP          3\n";
	ofile << "[ atoms ]\n";
	int atom_counter = 1;
	for(int j = 0; j < box.res.size(); j++)
	for(int i = 0; i < box.res[j].atoms.size(); i++)
	{
		ofile << atom_counter << "     " <<  box.res[j].atoms[i].type << "     " <<  j+1 << "     " <<  box.res[j].name << "     " <<  box.res[j].atoms[i].name;
  		ofile << "     " << atom_counter << "     " << box.res[j].atoms[i].charge << "     " << box.res[j].atoms[i].mass << "\n";
		atom_counter++;
	}
	ofile << "[ bonds ]\n";
	for(int i = 0; i < box.bonds.size(); i++)
	{
		ofile << box.bonds[i].i+1 << "     " <<  box.bonds[i].j+1 << "     " <<  box.bonds[i].function;
		if(box.bonds[i].has_params) ofile << "     " << box.bonds[i].b0 << "     " << box.bonds[i].kb;
		if(box.bonds[i].function == 1) ofile << "     ; Harmonic potential ";
		else if(box.bonds[i].function == 3) ofile << "     " << box.bonds[i].a << "     ; Morse potential ";
		else if(box.bonds[i].function == 6) ofile << "     ; Harmonic No Exclusions ";
		else ofile << "     ; Other Bond Type";
		ofile << box.atoms[box.bonds[i].i].name << "-" << box.atoms[box.bonds[i].j].name << "\n";
	}

	ofile << "[ pairs ]\n";
	std::cout << "exclusions.size == " << box.exclusions.size() << std::endl;
	for (int i=0; i< box.exclusions.size();i++)
	{
		ofile << box.exclusions[i].i << "  " << box.exclusions[i].j << "    " << box.exclusions[i].function << "\n";
	}

	ofile << "[ angles ]\n";
	for(int i = 0; i < box.angles.size(); i++)
	{
		ofile << box.angles[i].i+1 << "     " <<  box.angles[i].j+1 << "     " << box.angles[i].k+1 << "     " <<  box.angles[i].function;
		ofile << "     "  << box.angles[i].angle << "     " << box.angles[i].kb;
		ofile << "; " << box.atoms[box.angles[i].i].name << "-" << box.atoms[box.angles[i].j].name << "-" << box.atoms[box.angles[i].k].name;
		ofile << "\n";
	}	
	ofile << "[ dihedrals ]\n";
	for(int i = 0; i < box.dihedrals.size(); i++)
	{
		ofile << box.dihedrals[i].i+1 << "     " <<  box.dihedrals[i].j+1 << "     " <<  box.dihedrals[i].k+1 << "     " <<  box.dihedrals[i].l+1 << "     " <<  box.dihedrals[i].function;
		ofile << "     " << box.dihedrals[i].angle << "     " << box.dihedrals[i].kb << "     " << box.dihedrals[i].multiplicity;
		ofile << "; " << box.atoms[box.dihedrals[i].i].name << "-" << box.atoms[box.dihedrals[i].j].name << "-" << box.atoms[box.dihedrals[i].k].name << "-" << box.atoms[box.dihedrals[i].l].name << "\n";
	}
	ofile.close();	
		
	
	return;
}



void get_com_local(std::vector<Atom> atoms, float& x, float& y, float& z, float xref, float yref, float zref)
{
	float x_tot = 0.0f;
	float y_tot = 0.0f;
	float z_tot = 0.0f;
	int counter = 0;
	for(int i = 0; i < atoms.size(); i++)
	{
		if((atoms[i].pos[0] - xref)*(atoms[i].pos[0] - xref) + (atoms[i].pos[1] - yref)*(atoms[i].pos[1] - yref) + (atoms[i].pos[2] - zref)*(atoms[i].pos[2] - zref) < 1.0f)
		{
			x_tot += atoms[i].pos[0];
			y_tot += atoms[i].pos[1];
			z_tot += atoms[i].pos[2];
			counter++;
		}
	}
	
	x = x_tot * 1.0f/(float)counter;
	y = y_tot * 1.0f/(float)counter;
	z = z_tot * 1.0f/(float)counter;
	return;
}
static inline float get_dist(float a, float b, float c, float x, float y, float z)
{
	return sqrt((a-x)*(a-x) + (b-y)*(b-y) + (c-z)*(c-z));
}
void get_3_nearest_golds(std::vector<Atom>& atoms, Atom sulfur, int& i1, int& i2, int& i3, float& min1, float& min2, float& min3) 
{
	i1 = 0; i2 = 0; i3 = 0;
	float x = sulfur.pos[0];
	float y = sulfur.pos[1];
	float z = sulfur.pos[2];
	int ih1=0, ih2=0;
	float mh1=1e10, mh2=1e10;
	min1 = 1e10; min2 = 1e10; min3 = 1e10;
	for(int i = 0; i < atoms.size(); i++) //should give nearest 3 gold atoms to a test point
	{
		float ax = atoms[i].pos[0]; float ay = atoms[i].pos[1]; float az = atoms[i].pos[2];
		float r = get_dist(ax, ay, az, x, y, z);
		ih1 = i; mh1 = r; //load params into mh1 and ih1
		if( mh1 < min1 )
		{
			ih2 = i1; mh2 = min1; //place old values into holder
			min1 = mh1; i1 = ih1; //replace with new vals
			ih1 = ih2; mh1 = mh2; //put old values into mh1 to ensure that it gets tested against other nearby golds before being removed from the list
		}
		if( mh1 < min2 )
		{
			ih2 = i2; mh2 = min2;
			min2 = mh1; i2 = ih1;
			ih1 = ih2; mh1 = mh2;
		}
		if( mh1 < min3 )
		{
			ih2 = i3; mh2 = min3;
			min3 = mh1; i3 = ih1;
			ih1 = ih2; mh1 = mh2;
		}
	}
	
	
	return;
}

void get_nearby_gold_indexes(std::vector<Atom>& atoms, Ligand lig, std::vector<int>& indexes) 
{
	indexes.clear();
	float x = lig.atoms[0].pos[0];
	float y = lig.atoms[0].pos[1];
	float z = lig.atoms[0].pos[2];
	for(int i = 0; i < atoms.size(); i++) //should give nearest 3 gold atoms to a test point
	{
		float ax = atoms[i].pos[0]; float ay = atoms[i].pos[1]; float az = atoms[i].pos[2];
		float r = get_dist(ax, ay, az, x, y, z);
		if(r < 1.0f)
		{
			indexes.push_back(i);
		}
	}
	
	
	return;
}

void build_nanoparticle(std::vector<Atom> gnp_atoms, std::vector<Atom> sulfurs, std::string ligand_top, std::string ligand_gro, std::string output_name, bool nm_flag)
{
	
	Box box;
	Residue ligand_res;
	Ligand ligand_template;
	ligand_template = load_ligand(ligand_top, ligand_gro);
	//ligand_template = make_c12_ligand();
	Residue gnp_res;
	gnp_res.name = "NP";
	gnp_res.index = 0;
	gnp_res.atoms = gnp_atoms;
	if(!nm_flag) generate_bonds(gnp_res, 0.75);
	else generate_bonds(gnp_res, 1.0);
	box.res.push_back(gnp_res);
	
	int atom_counter= gnp_res.atoms.size();
	int ligand_counter = 1;
	
	for(int i = 0; i < sulfurs.size(); i++)
	{
		Ligand ligand = ligand_template;
		ligand.index = ligand_counter;
		ligand_counter++;
		float xcom,ycom,zcom;
		get_com_local(gnp_atoms, xcom, ycom, zcom, sulfurs[i].pos[0], sulfurs[i].pos[1], sulfurs[i].pos[2]);
		Eigen::Vector3d vaxis;
		vaxis << sulfurs[i].pos[0] - xcom, sulfurs[i].pos[1] - ycom, sulfurs[i].pos[2] - zcom;
		vaxis.normalize();
		ligand.rotate_to(vaxis);
		//std::cout << ligand.atoms[0].pos << std::endl;
		for(int j = 0; j < ligand.atoms.size(); j++) //create a ligand for every sulfur
		{
			ligand.atoms[j].pos[0] += sulfurs[i].pos[0];
			ligand.atoms[j].pos[1] += sulfurs[i].pos[1];
			ligand.atoms[j].pos[2] += sulfurs[i].pos[2];
		}
		box.res.push_back(lig2res(ligand));
	}
	atom_counter = 0;
	for(int i = 0; i < box.res.size(); i++)
	{
		box.res[i].index = i;
		box.res[i].offset = atom_counter;
		for(int j = 0; j < box.res[i].atoms.size(); j++)
		{
			box.res[i].atoms[j].index = atom_counter;
			atom_counter++;
		}
	}
	
	for(int i = 1; i < box.res.size(); i++)
	{
		int l_index = box.res[i].offset;
		int i1, i2, i3;
		float min1, min2, min3;
		get_3_nearest_golds(box.res[0].atoms, box.res[i].atoms[0], i1, i2, i3, min1, min2, min3);
		Bond b1;
		b1.function = 6; //harmonic, no exclusions
		b1.has_params = 1;
		b1.i = i1;
		b1.j = l_index;
		b1.b0 = min1;
		b1.kb = (14.7*14.7)*2*36.664;
		box.bonds.push_back(b1);
		b1.i = i2;
		b1.j = l_index;
		b1.b0 = min2;
		b1.kb = (14.7*14.7)*2*36.664;
		box.bonds.push_back(b1);
		b1.i = i3;
		b1.j = l_index;
		b1.b0 = min3;
		b1.kb = (14.7*14.7)*2*36.664;
		box.bonds.push_back(b1);
		ligand_counter++;	
	}
	for(int i = 0; i < box.res.size(); i++)
	{
		for(int j = 0; j < box.res[i].atoms.size(); j++)
		{
			box.atoms.push_back(box.res[i].atoms[j]);
		}
	}
	//combine bonds from various residue
	std::vector<Bond> bonds;
	std::vector<Pair> exclusions;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;
	for(int i = 0; i < box.res.size(); i++)
	{
		for(int j = 0; j < box.res[i].bonds.size(); j++)
		{
			Bond bcopy = box.res[i].bonds[j];
			bcopy.i += box.res[i].offset;
			bcopy.j += box.res[i].offset;
			bonds.push_back(bcopy);
		}
		for(int j = 0; j < box.res[i].exclusions.size(); j++)
		{
			Pair bcopy = box.res[i].exclusions[j];
			bcopy.i += box.res[i].offset;
			bcopy.j += box.res[i].offset;
			exclusions.push_back(bcopy);
		}
		for(int j = 0; j < box.res[i].angles.size(); j++)
		{
			Angle bcopy = box.res[i].angles[j];
			bcopy.i += box.res[i].offset;
			bcopy.j += box.res[i].offset;
			bcopy.k += box.res[i].offset;
			angles.push_back(bcopy);
		}
		for(int j = 0; j < box.res[i].dihedrals.size(); j++)
		{
			Dihedral bcopy = box.res[i].dihedrals[j];
			bcopy.i += box.res[i].offset;
			bcopy.j += box.res[i].offset;
			bcopy.k += box.res[i].offset;
			bcopy.l += box.res[i].offset;
			dihedrals.push_back(bcopy);
		}
		
	}
	box.bonds.insert( box.bonds.end(), bonds.begin(), bonds.end() );
	box.exclusions.insert( box.exclusions.end(), exclusions.begin(), exclusions.end() );
	box.angles.insert( box.angles.end(), angles.begin(), angles.end() );
	box.dihedrals.insert( box.dihedrals.end(), dihedrals.begin(), dihedrals.end() );
	write_itp(output_name, box);
	dump_atoms_gro(output_name, box);	
	return;
}

int main(int argc, char **argv)
{
	// the input name
	std::string fname = argv[1];
	InputParser parser;

	ParameterPack pack_;
	parser.ParseFile(fname, pack_);

	auto ligandPack = pack_.findParamPack("ligand", ParameterPack::KeyType::Required);
	auto NPPack		= pack_.findParamPack("Nanoparticle", ParameterPack::KeyType::Optional);
	auto ShapePack  = pack_.findParamPack("Shape", ParameterPack::KeyType::Optional);
	auto OutPack 	= pack_.findParamPack("output", ParameterPack::KeyType::Optional);

	//for ligand
	std::string ligand_top = "";
	std::string ligand_gro = "";
	ligandPack->ReadString("topology", ParameterPack::KeyType::Required, ligand_top);
	ligandPack->ReadNumber("gro", ParameterPack::KeyType::Required, ligand_gro);

	//for manually-specified nanoparticle
	std::string np_filename = "";
	int np_assigned = 0;
	if (NPPack != nullptr)
	{
		np_assigned = 1;
		NPPack -> ReadString("name", ParameterPack::KeyType::Required, np_filename);
	}

	//for shape
	std::string np_shape;
	float params[2];
	//for auto-generated nanoparticle
	bool generate_np = 0;
	bool np_file_specified = 0;

	if (ShapePack != nullptr)
	{
		generate_np = 1;

		ShapePack -> ReadString("type", ParameterPack::KeyType::Required, np_shape);

		if (np_shape =="sphere")
		{
			float radius;
			ShapePack -> ReadNumber("radius", ParameterPack::KeyType::Required, radius);

			params[0] = radius;
		}

		if (np_shape=="cylinder")
		{
			float radius;
			float length;
			ShapePack -> ReadNumber("radius", ParameterPack::KeyType::Required, radius);
			ShapePack -> ReadNumber("length", ParameterPack::KeyType::Required, length);

			params[0] = radius;
			params[1] = length;
		}
	}

	// read output things
	std::string output_name = "out";
	float e = 50.0; //kJ/mol
	if (OutPack != nullptr)
	{
		OutPack -> ReadNumber("energy", ParameterPack::KeyType::Optional,e);
		OutPack -> ReadString("name", ParameterPack::KeyType::Optional, output_name);
	}
	
	std::vector<Atom> atoms;
	std::vector<Atom> lig_ats;
	std::vector<Atom> sulfurs;
	bool skip_grid = 0;
	if(generate_np)
	{
		gen_lattice_fcc(4*params[0],4*params[0],4*params[0],atoms);
		if(np_shape == "sphere"){
			cut_sphere(params[0], 0.4065*round(2.0f*params[0]/0.4065), 0.4065*round(2.0f*params[0]/0.4065), 0.4065*round(2.0f*params[0]/0.4065), atoms);
		}
		else if(np_shape == "cylinder"){
			cut_fcc_above(4*params[0] + params[1]/2.0f, atoms);
			cut_fcc_below(4*params[0] - params[1]/2.0f, atoms);
			float p1[3] = {0,0,0};
			float p2[3] = {1,1,1};
			float p3[3] = {0,1,0};
			cut_cylinder(params[0], p1, p2, atoms);
			rotate_atoms_about_com(atoms, p2, p3);
		}
		else return 0;
	}
	else if(np_assigned == 1){
		read_gro_np(np_filename, atoms);
	}
	else if(np_assigned == 2){
		read_gro_np2(np_filename, atoms, sulfurs); //puts gold atom types in atoms, puts sulfur atom types in sulfur, skips the rest, potential grid skipped
		skip_grid = 1;
	}
	else return 0;
	std::cout << "Number of gold atoms = " << atoms.size() << std::endl;
	if(!skip_grid) generate_potential_grid(atoms, sulfurs, e);
	build_nanoparticle(atoms,sulfurs, ligand_top, ligand_gro, output_name, 1);

	return 0;
}
