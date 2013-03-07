
extern "C" int mainClassify(int argc, char *argv[]);

extern "C" int mainRegPath(int argc, char *argv[]);

extern "C" int mainTain(int argc, char *argv[]);

int main(int argc, char* argv[])
{
	argc = 6;
	char** argvloc = new char*[9];

	argvloc[0]= "";
	argvloc[1]= "-s";
	argvloc[2]= "D:\\Packages\\Ell1_Logistic\\examples\\exs_internetad_X";
	argvloc[3]= "D:\\Packages\\Ell1_Logistic\\examples\\exs_internetad_b";
	argvloc[4]= "0.01";
	argvloc[5]= "D:\\Packages\\Ell1_Logistic\\examples\\model_internetad";

	int resultTain = mainTain( argc, argvloc);
	
	argvloc[0]= "";
	argvloc[1]= "-t";
	argvloc[2]= "D:\\Packages\\Ell1_Logistic\\examples\\exs_internetad_b";
	argvloc[3]= "D:\\Packages\\Ell1_Logistic\\examples\\model_internetad";
	argvloc[4]= "D:\\Packages\\Ell1_Logistic\\examples\\exs_internetad_X";
	argvloc[5]= "D:\\Packages\\Ell1_Logistic\\examples\\result_internetad";
	//-t exs_internetad_b model_internetad exs_internetad_X result_internetad
	int resultClassify = mainClassify(argc, argvloc);

	argvloc[0]= "";
	argvloc[1]= "-r";
	argvloc[2]= "-s";
	argvloc[3]= "D:\\Packages\\Ell1_Logistic\\examples\\exs_internetad_X";
	argvloc[4]= "D:\\Packages\\Ell1_Logistic\\examples\\exs_internetad_b";
	argvloc[5]= "0.00001";//Was "0.00001"
	argvloc[6]= "1000";  //Was 100
	argvloc[7]= "D:\\Packages\\Ell1_Logistic\\examples\\path_internetad";
	argvloc[8]= "-v";
	argvloc[9]= "3";
	//-r -s exs_internetad_X exs_internetad_b 0.0001 100 path_internetad -v 3
	int resultRegPath = mainRegPath( argc, argvloc);

	
	return 0;
}

