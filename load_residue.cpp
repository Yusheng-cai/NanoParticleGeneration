#include "basic_structs.h"
// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}


static inline bool check_comment(std::string line)
{
	if(line.length() == 0) return true;
	std::stringstream ss(line);
	std::string word;
	ss >> word;

	if(*word.c_str() == ';' || *word.c_str() == '#' || *word.c_str() == '%' || line.length() == 0) return true;
	return false;
}

static std::string load_atoms(std::ifstream& ifile, std::vector<Atom>& atoms)
{
	std::string line;
	Atom a1;
	while(getline(ifile, line))
	{
		if(check_comment(line)) continue;
		std::stringstream ss(line);
		std::string test;
		ss >> test;

		if (test == "[")
		{
			std::cout << "line = " << line << std::endl;
			break;
		}

		ss = std::stringstream(line);
		ss >> a1.index >> a1.type >> a1.resnr >> a1.resname >> a1.name >> a1.chargegrp >> a1.charge >> a1.mass;
		a1.index--;
		a1.resnr--;
		atoms.push_back(a1);
	}
	return line;
}

static std::string load_bonds(std::ifstream& ifile, std::vector<Bond>& bonds)
{
	std::string line;
	Bond b1;
	while(getline(ifile, line))
	{
		if(check_comment(line)) continue;
		if(line.at(0) == '[') break;
		std::stringstream ss(line);
		ss >> b1.i >> b1.j >> b1.function >> b1.b0 >> b1.kb;
		b1.i--;
		b1.j--;
		bonds.push_back(b1);
	}
	return line;	
}

static std::string load_pairs(std::ifstream& ifile, std::vector<Pair>& pairs)
{
	std::string line;
	Pair p1;
	while(getline(ifile, line))
	{
		if(check_comment(line)) continue;
		if(line.at(0) == '[') break;
		std::stringstream ss(line);
		ss >> p1.i >> p1.j >> p1.function;
		p1.i--;
		p1.j--;
		pairs.push_back(p1);
	}
	return line;		
}

static std::string load_angles(std::ifstream& ifile, std::vector<Angle>& angles)
{
	std::string line;
	Angle a1;
	while(getline(ifile, line))
	{
		if(check_comment(line)) continue;
		if(line.at(0) == '[') break;
		std::stringstream ss(line);
		ss >> a1.i >> a1.j >> a1.k >> a1.function >> a1.angle >> a1.kb;
		a1.i--;
		a1.j--;
		a1.k--;
		angles.push_back(a1);
	}
	return line;		
}

static std::string load_dihedrals(std::ifstream& ifile, std::vector<Dihedral>& dihedrals)
{
	std::string line;
	Dihedral d1;
	while(getline(ifile, line))
	{
		if(check_comment(line)) continue;
		if(line.at(0) == '[') break;
		std::stringstream ss(line);
		ss >> d1.i >> d1.j >> d1.k >> d1.l >> d1.function >> d1.angle >> d1.kb >> d1.multiplicity;
		d1.i--;
		d1.j--;
		d1.k--;
		d1.l--;
		dihedrals.push_back(d1);
	}
	return line;			
}
static int read_top(std::string topology_filename, Ligand& res)
{
	std::ifstream ifile(topology_filename);
	if(!ifile.is_open())
	{
		std::cout << "Failed to open topology file for ligand." << std::endl;
		return -1;
	}

	std::string line;
	Atom a_holder;

	while(getline(ifile, line))
	{
		if(check_comment(line)) continue;
		if(line.find("[ moleculetype ]") != std::string::npos){
			while(true){
				getline(ifile, line);
				if(check_comment(line)) continue;
				std::stringstream ss(line);
				ss >> res.name >> res.nrexcl;
				break;
			}
			std::stringstream ss;
			while(getline(ifile, line)){
				if(line.find("[ atoms ]") != std::string::npos) line = load_atoms(ifile, res.atoms);
				if(line.find("[ bonds ]") != std::string::npos) line = load_bonds(ifile, res.bonds);
				if(line.find("[ pairs ]") != std::string::npos) line = load_pairs(ifile, res.exclusions);
				if(line.find("[ angles ]") != std::string::npos) line = load_angles(ifile, res.angles);
				if(line.find("[ dihedrals ]") != std::string::npos) line = load_dihedrals(ifile, res.dihedrals);
			}
			break;
		}
	}
	return 1;
}

static int read_gro(std::string filename, Ligand& res)
{
	std::ifstream ifile(filename);
	if(!ifile.is_open())
	{
		std::cout << "Failed to open gro file for ligand." << std::endl;
		return -1;
	}
	std::string line;
	getline(ifile, line);
	getline(ifile, line);
	int num_atoms = std::stoi(line);
	if(num_atoms != res.atoms.size())
	{
		std::cout << "Disagreement in atom count between gro file and top file." << std::endl;
		return -1;
	}
	int atom_counter = 0;
	while(atom_counter < num_atoms)
	{
		getline(ifile, line);
		if(check_comment(line)) continue;
		res.atoms[atom_counter].pos[0] = std::stof(line.substr(20, 8));
		res.atoms[atom_counter].pos[1] = std::stof(line.substr(28, 8));
		res.atoms[atom_counter].pos[2] = std::stof(line.substr(36, 8));
		res.atoms[atom_counter].vx = 0.0f;
		res.atoms[atom_counter].vy = 0.0f;
		res.atoms[atom_counter].vz = 0.0f;
		atom_counter++;
	}
	ifile.close();
	return 1;
}

int read_gro_np(std::string filename, std::vector<Atom>& atoms)
{
	std::ifstream ifile(filename);
	if(!ifile.is_open())
	{
		std::cout << "Failed to open gro file for nanoparticle." << std::endl;
		return -1;
	}
	std::string line;
	getline(ifile, line);
	getline(ifile, line);
	int num_atoms = std::stoi(line);
	int atom_counter = 0;
	while(atom_counter < num_atoms)
	{
		getline(ifile, line);
		if(check_comment(line)) continue;
		Atom a1;
		a1.pos[0] = std::stof(line.substr(20, 8));
		a1.pos[1] = std::stof(line.substr(28, 8));
		a1.pos[2] = std::stof(line.substr(36, 8));
		a1.vx = 0.0f;
		a1.vy = 0.0f;
		a1.vz = 0.0f;
		a1.name = line.substr(10,5);
		trim(a1.name);
		if(a1.name == "AU" || a1.name == "AUS" || a1.name == "AUL"){
			a1.name = "AU";
			a1.type = "AU";
			a1.mass = 196.96;
			a1.charge = 0.0f;
			atoms.push_back(a1);
		}
		atom_counter++;
	}
	ifile.close();
	return 1;
}

int read_gro_np2(std::string filename, std::vector<Atom>& atoms, std::vector<Atom>& sulfurs)
{
	std::ifstream ifile(filename);
	if(!ifile.is_open())
	{
		std::cout << "Failed to open gro file for nanoparticle." << std::endl;
		return -1;
	}
	std::string line;
	getline(ifile, line);
	getline(ifile, line);
	int num_atoms = std::stoi(line);
	int atom_counter = 0;
	while(atom_counter < num_atoms)
	{
		getline(ifile, line);
		if(check_comment(line)) continue;
		Atom a1;
		a1.pos[0] = std::stof(line.substr(20, 8));
		a1.pos[1] = std::stof(line.substr(28, 8));
		a1.pos[2] = std::stof(line.substr(36, 8));
		a1.vx = 0.0f;
		a1.vy = 0.0f;
		a1.vz = 0.0f;
		a1.name = line.substr(10,5);
		trim(a1.name);
		if(a1.name == "AU" || a1.name == "AUS" || a1.name == "AUL"){
			a1.name = "AU";
			a1.type = "AU";
			a1.mass = 196.96;
			a1.charge = 0.0f;
			atoms.push_back(a1);
		}
		if(a1.name.at(0) == 'S'){
			a1.name = "S";
			a1.type = "S";
			a1.mass = 32.06f;
			a1.charge = 0.0f;
			sulfurs.push_back(a1);
		}
		atom_counter++;		
	}
	ifile.close();
	return 1;
}
Ligand load_ligand(std::string topology_filename, std::string gro_filename)
{
	Ligand res;
	read_top(topology_filename, res);
	read_gro(gro_filename, res);
	res.shift(-res.atoms[0].pos[0], -res.atoms[0].pos[1], -res.atoms[0].pos[2]); //move sulfur to origin
	res.set_axis();
	return res;
}



Ligand make_c18_ligand()
{
	Ligand ligand;
	ligand.add(Atom("S",      -0.64612473016,   -0.02529243595,   -0.00321153372));                 
	ligand.add(Atom("U1",     -0.49913025094,    0.08085896397,    0.00668245337));                 
	ligand.add(Atom("U1",     -0.37169844144,   -0.00257302108,   -0.00269728189));                 
	ligand.add(Atom("U1",     -0.24727953449,    0.08634214890,    0.00501698876));                 
	ligand.add(Atom("U1",     -0.11913176614,    0.00334792059,   -0.00409570121));                 
	ligand.add(Atom("U1",      0.00546008745,    0.09177103054,    0.00363366598));                 
	ligand.add(Atom("U1",      0.13317595396,    0.00810043607,   -0.00523094309));                 
	ligand.add(Atom("U1",      0.25834989630,    0.09568040222,    0.00261795984));                 
	ligand.add(Atom("U1",      0.38540588349,    0.01096524392,   -0.00592551027));                 
	ligand.add(Atom("U1",      0.51135152829,    0.09739662094,    0.00213333594));                 
	ligand.add(Atom("U1",      0.63759899948,    0.01142288149,   -0.00600091273));                 
	ligand.add(Atom("U1",      0.76437583739,    0.09658353772,    0.00233896998));                 
	ligand.add(Atom("U1",      0.88980320833,    0.00935513405,   -0.00530982402));                 
	ligand.add(Atom("U1",      1.01732143984,    0.09335362602,    0.00333544923));                 
	ligand.add(Atom("U1",      1.14206656301,    0.00509751892,   -0.00381910231));                 
	ligand.add(Atom("U1",      1.27011917525,    0.08824783982,    0.00507588266));                 
	ligand.add(Atom("U1",      1.39443062747,   -0.00064651261,   -0.00169067208));                 
	ligand.add(Atom("U1",      1.52265785981,    0.08202065243,    0.00731624881));                 
	ligand.add(Atom("CT",      1.64630561734,   -0.00622336359,    0.00081270476));                 
	ligand.add(Atom("HC",      1.64765508513,   -0.07801504050,    0.08344505088));                 
	ligand.add(Atom("HC",      1.73687750974,    0.05486531482,    0.00741178214));                 
	ligand.add(Atom("HC",      1.64975751276,   -0.06205376281,   -0.09328691495));  
	ligand.shift(0.64612473016, 0.02529243595, 0.00321153372); //move sulfur to origin
	ligand.set_axis();

	ligand.add_bond(0,1); //S-U
	ligand.add_bond(1,2); //U-U
	ligand.add_bond(2,3); //U-U
	ligand.add_bond(3,4); //U-U
	ligand.add_bond(4,5); //U-U
	ligand.add_bond(5,6); //U-U
	ligand.add_bond(6,7); //U-U
	ligand.add_bond(7,8); //U-U
	ligand.add_bond(8,9);  //U-U	
	ligand.add_bond(9,10); //U-U
	ligand.add_bond(10,11); //U-U
	ligand.add_bond(11,12); //U-U
	ligand.add_bond(12,13); //U-U
	ligand.add_bond(13,14); //U-U
	ligand.add_bond(14,15);  //U-U		
	ligand.add_bond(15,16); //U-C
	ligand.add_bond(16,17); //U-C
	ligand.add_bond(17,18); //U-C
	ligand.add_bond(18,19); //H-C
	ligand.add_bond(18,20); //H-C
	ligand.add_bond(18,21); //H-C

	ligand.add_angle(0,1,2); //S-U-U
	ligand.add_angle(1,2,3); //U-U-U
	ligand.add_angle(2,3,4); //U-U-U
	ligand.add_angle(3,4,5); //U-U-U
	ligand.add_angle(4,5,6); //U-U-U
	ligand.add_angle(5,6,7); //U-U-U
	ligand.add_angle(6,7,8); //U-U-U
	ligand.add_angle(7,8,9); //U-U-U
	ligand.add_angle(8,9,10);  //U-U-U
	ligand.add_angle(9,10,11);  //U-U-U	
	ligand.add_angle(10,11,12);  //U-U-U
	ligand.add_angle(11,12,13); //U-U-U
	ligand.add_angle(12,13,14); //U-U-U
	ligand.add_angle(13,14,15); //U-U-U
	ligand.add_angle(14,15,16); //U-U-U
	ligand.add_angle(15,16,17);  //U-U-U
	ligand.add_angle(16,17,18);  //U-U-CT	
	ligand.add_angle(17,18,19); //U-CT-HC
	ligand.add_angle(17,18,20); //U-CT-HC
	ligand.add_angle(17,18,21); //U-CT-HC
	ligand.add_angle(19,18,20); //H-CT-HC
	ligand.add_angle(20,18,21); //H-CT-HC
	ligand.add_angle(19,18,21); //H-CT-HC
	
	ligand.add_dihedral(0,1,2,3); //S-U-U-U
	ligand.add_dihedral(1,2,3,4); //U-U-U-U
	ligand.add_dihedral(2,3,4,5); //U-U-U-U
	ligand.add_dihedral(3,4,5,6); //U-U-U-U
	ligand.add_dihedral(4,5,6,7); //U-U-U-U
	ligand.add_dihedral(5,6,7,8); //U-U-U-U
	ligand.add_dihedral(6,7,8,9); //U-U-U-U
	ligand.add_dihedral(7,8,9,10); //U-U-U-U
	ligand.add_dihedral(8,9,10,11);  //U-U-U-U
	ligand.add_dihedral(9,10,11,12); //U-U-U-U
	ligand.add_dihedral(10,11,12,13); //U-U-U-U
	ligand.add_dihedral(11,12,13,14); //U-U-U-U
	ligand.add_dihedral(12,13,14,15); //U-U-U-U
	ligand.add_dihedral(13,14,15,16); //U-U-U-U
	ligand.add_dihedral(14,15,16,17);  //U-U-U-U
	ligand.add_dihedral(15,16,17,18); //U-U-U-C
	ligand.add_dihedral(16,17,18,19); //U-U-C-H
	ligand.add_dihedral(16,17,18,20); //U-U-C-H
	ligand.add_dihedral(16,17,18,21); //U-U-C-H

	return ligand;
}
Ligand make_c12_ligand()
{
	Ligand ligand;
	ligand.add(Atom("S", -0.5998913, 0.05802054394, -0.00192137829));     //0
	ligand.add(Atom("U1", -0.44193055428, 0.14716011754,  0.00664366240));//1
	ligand.add(Atom("U1", -0.32462489424, 0.04998231971, -0.00264876992));//2
	ligand.add(Atom("U1", -0.19109354754, 0.12457941088,  0.00413923427));//3
	ligand.add(Atom("U1", -0.07285945951, 0.02793643419, -0.00472292315));//4          
	ligand.add(Atom("U1", 0.06057362590,  0.10236756561,  0.00249568077));//5
	ligand.add(Atom("U1", 0.17870895206,  0.00555226806, -0.00566466423));//6
	ligand.add(Atom("U1", 0.31221880774,  0.07977657684,  0.00219886480));//7
	ligand.add(Atom("U1", 0.43024142906, -0.01724584480, -0.00514976375));//8
	ligand.add(Atom("U1", 0.56382199117,  0.05678501512,  0.00330863274));//9
	ligand.add(Atom("U1", 0.68175560800, -0.04039103654, -0.00340438075));//10        	                 
	ligand.add(Atom("U1", 0.81527410816,  0.03348145143,  0.00536530963));//11  
	ligand.add(Atom("CT", 0.93273562980, -0.06284152714, -0.00098524383));//12  
	ligand.add(Atom("HC", 0.92932314534, -0.13437953424,  0.08180598331));//13  
	ligand.add(Atom("HC", 1.02718675611, -0.00791916241,  0.00545945660));//14  
	ligand.add(Atom("HC", 0.93241177585, -0.11897058716, -0.09496721796));//15
	ligand.shift(0.5998913, -0.05802054394, 0.00192137829); //move sulfur to origin
	ligand.set_axis();

	ligand.add_bond(0,1); //S-U
	ligand.add_bond(1,2); //U-U
	ligand.add_bond(2,3); //U-U
	ligand.add_bond(3,4); //U-U
	ligand.add_bond(4,5); //U-U
	ligand.add_bond(5,6); //U-U
	ligand.add_bond(6,7); //U-U
	ligand.add_bond(7,8); //U-U
	ligand.add_bond(8,9);  //U-U	
	ligand.add_bond(9,10); //U-C
	ligand.add_bond(10,11); //U-C
	ligand.add_bond(11,12); //U-C
	ligand.add_bond(12,13); //H-C
	ligand.add_bond(12,14); //H-C
	ligand.add_bond(12,15); //H-C

	ligand.add_angle(0,1,2); //S-U-U
	ligand.add_angle(1,2,3); //U-U-U
	ligand.add_angle(2,3,4); //U-U-U
	ligand.add_angle(3,4,5); //U-U-U
	ligand.add_angle(4,5,6); //U-U-U
	ligand.add_angle(5,6,7); //U-U-U
	ligand.add_angle(6,7,8); //U-U-U
	ligand.add_angle(7,8,9); //U-U-U
	ligand.add_angle(8,9,10);  //U-U-U
	ligand.add_angle(9,10,11);  //U-U-U	
	ligand.add_angle(10,11,12);  //U-U-U		
	ligand.add_angle(11,12,13); //U-C-H
	ligand.add_angle(11,12,14); //U-C-H
	ligand.add_angle(11,12,15); //U-C-H
	ligand.add_angle(13,12,14); //H-C-H
	ligand.add_angle(14,12,15); //H-C-H
	ligand.add_angle(13,12,15); //H-C-H
	
	ligand.add_dihedral(0,1,2,3); //S-U-U-U
	ligand.add_dihedral(1,2,3,4); //U-U-U-U
	ligand.add_dihedral(2,3,4,5); //U-U-U-U
	ligand.add_dihedral(3,4,5,6); //U-U-U-U
	ligand.add_dihedral(4,5,6,7); //U-U-U-U
	ligand.add_dihedral(5,6,7,8); //U-U-U-U
	ligand.add_dihedral(6,7,8,9); //U-U-U-U
	ligand.add_dihedral(7,8,9,10); //U-U-U-U
	ligand.add_dihedral(8,9,10,11);  //U-U-U-U
	ligand.add_dihedral(9,10,11,12); //U-U-U-C
	ligand.add_dihedral(10,11,12,13); //U-U-C-H
	ligand.add_dihedral(10,11,12,14); //U-U-C-H
	ligand.add_dihedral(10,11,12,15); //U-U-C-H

	return ligand;
}