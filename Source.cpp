#include <iostream>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <ctime>
using namespace std;


double* read(char x[10])
{
	double conf[6]; 
	ifstream fin(x);
	fin >> conf[1] >> conf[2] >> conf[3] >> conf[4] >> conf [5];
	fin.close();
	return conf;
}


bool distcheck(double x1, double x2, double y1, double y2, double z1, double z2, double molrad, double length)
{
	double r;
	double delta_x = abs(x2-x1), delta_y = abs(y2 - y1), delta_z = abs(z2 - z1);
	if (delta_x > length/2)
	{
		delta_x = length - delta_x;
	}
	if (delta_y > length/2)
	{
		delta_y = length - delta_y;
	}
	if (delta_z > length/2)
	{
		delta_z = length - delta_z;
	}
	r = sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
	if (r > 2 * molrad)
	{
		return (true);
	}
	else
	{
		return (false);
	}
}

double jump(double coor,double dcoor, double length)
{
	coor = coor + (rand() * 2.0 / RAND_MAX) * dcoor - dcoor;
	if (coor >= length)
	{
		coor = coor - length;
	}
	else
	{
		if (coor < 0)
		{
			coor = coor + length;
		}
	}
	return (coor);
}


int main()
{
	srand(static_cast<unsigned int>(time(0)));
	char way[10]; double partnum, steps, length, dcoor, molrad; int badsteps=0; int shots = 0; double sigma = 3.4 * pow(10, -10);
	cin >> way;
	partnum = read(way)[1];
	steps = read(way)[2];
	length = read(way)[3]* sigma;
	dcoor = read(way)[4]* sigma;
	molrad = read(way)[5]* sigma;
	double* Xcoor = new double[partnum + 1];
	double* Ycoor = new double[partnum + 1];
	double* Zcoor = new double[partnum + 1];
	int n = ceil(pow(partnum, 1.0 / 3)); // количество узлов в сетке на грань куба
	double griddist = length / (n);
	// cout << griddist << endl;
	
	for (int i = 1; i <= partnum; i++)
	{
		Xcoor[i] = griddist * ((i % n));
		Ycoor[i] = griddist * (((int) (i / n)) % n);
		Zcoor[i] = griddist * ((int) (i / (n*n)) % n);
	}
	ofstream fout;
	fout.open("1.XMOL");
	fout << partnum << endl;
	fout << "STEP=" << 0 << endl;
	for (int i = 1; i <= partnum; i++)
	{
		fout << "Ar " << Xcoor[i] << " " << Ycoor[i] << " " << Zcoor[i] << endl;
	}
	
	for (int stepcount = 1; stepcount <= steps; stepcount++)
	{
		int r; 
		r = rand() % int (partnum) + 1; //выбрали частицу
		

		double Xcoorrem, Ycoorrem, Zcoorrem;
		Xcoorrem = Xcoor[r]; Ycoorrem = Ycoor[r]; Zcoorrem = Zcoor[r]; //запомнили положение частицы
		
		Xcoor[r] = jump(Xcoor[r], dcoor, length);
		Ycoor[r] = jump(Ycoor[r], dcoor, length); //прыжок частицей
		Zcoor[r] = jump(Zcoor[r], dcoor, length);
		
		bool t=true; //флаг принятия шага
		
		for (int k = 1; k <= partnum; k++) //проверка удалённости других частиц
		{	
			//cout << "step=" << stepcount << " part1=" << r << " part2=" << k << endl;
			if (k == r)
			{
				t = true; 
			}
			else
			{
				t = distcheck(Xcoor[k], Xcoor[r], Ycoor[k], Ycoor[r], Zcoor[k], Zcoor[r], molrad, length); 
			}
			//cout << t << endl;

			if (t == false)
			{
				//cout << "OOOOPS!" << endl << endl << endl; 
				break;
			}
		}

		if (t == false) //принятие или непринятие шага
		{
			badsteps = badsteps + 1;
			Xcoor[r] = Xcoorrem;
			Ycoor[r] = Ycoorrem;
			Zcoor[r] = Zcoorrem;
		}
		
		if ((stepcount % 125) == 0)
		{
			fout << partnum << endl;
			fout << "STEP=" << stepcount << endl;
			shots = shots + 1;
			for (int i = 1; i <= partnum; i++)
			{
				fout << "Ar " << Xcoor[i] << " " << Ycoor[i] << " " << Zcoor[i] << endl;
			}
		}
	}
	fout.close();
	fout.open("c1.txt");
	fout << partnum<<" "<< steps<<" "<< length<<" "<< dcoor<<" "<< molrad<<" "<<shots<< endl;
	fout.close();
	double t = badsteps / steps;
	cout << "bad = " << t << endl;
	cout << endl << partnum << " " << steps << " " << length << " " << dcoor << " " << molrad << endl;
	system("pause");
	return 0;
}
