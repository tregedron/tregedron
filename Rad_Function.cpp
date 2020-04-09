#include <iostream>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <string>

using namespace std;

const double pi = 3.1415926535897932;

double* read(char x[10]) //прочтение файла конфигурации
{
	double conf[7];
	ifstream fin(x);
	fin >> conf[1] >> conf[2] >> conf[3] >> conf[4] >> conf[5] >> conf[6];
	fin.close();
	return conf;
}

double dist(double x1, double x2, double y1, double y2, double z1, double z2, double length) //дистанция между частицами
{
	double r;
	double delta_x = abs(x2 - x1), delta_y = abs(y2 - y1), delta_z = abs(z2 - z1);
	if (delta_x > length / 2)
	{
		delta_x = length - delta_x;
	}
	if (delta_y > length / 2)
	{
		delta_y = length - delta_y;
	}
	if (delta_z > length / 2)
	{
		delta_z = length - delta_z;
	}
	r = sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
	return r;
}

double raspred(double n, double i, double delta, double norm) //функция для расчёта по локальной плотности функции распределения
{
	double gr = 0; 
	double R = (i + 0.5) * delta;
	gr = n / (4 * pi * pow(R, 2) * delta) / norm;
	return gr;
}

int main()
{
	char config[10], coordinatefile[10]; double partnum, length, shots, sloi, delta, pairdensity; int fshots;
	cin >> config;
	sloi = 120;
	cin >> fshots;
	partnum = read(config)[1];
	length = read(config)[3];
	shots = read(config)[6];
	delta = length / (2 * sloi);
	pairdensity = (partnum * (partnum - 1) / 2)/(pow(length,3));
	cout << endl << partnum << " " << shots << " " << length << endl;
	cin >> coordinatefile;
	double* Xcoor = new double[partnum];
	double* Ycoor = new double[partnum];
	double* Zcoor = new double[partnum];
	double* g = new double[sloi];
	double* gmass = new double[sloi];
	string line = "";
	string lineX = ""; string lineY = ""; string lineZ = "";
	ifstream fin(coordinatefile);
	ofstream fout("gABS.txt");
	for (int i = 0; i < sloi; i++)
	{
		g[i] = 0; gmass[i] = 0;
	}
	for (int i = 0; i < shots; i++) //прочтение файла положений
	{
		for (int k = 1; k <= partnum + 2; k++)
		{
			getline(fin, line);
			if (k >= 3)
			{
				int size = line.size();
				int pos = line.find(" ") + 1;
				lineX = line.substr(pos, size);
				size = lineX.size();
				pos = lineX.find(" ") + 1;
				lineY = lineX.substr(pos, size);
				size = lineY.size();
				pos = lineY.find(" ") + 1;
				lineZ = lineY.substr(pos, size);
				pos = lineX.find(" ") + 1;
				lineX = lineX.substr(0, pos);
				pos = lineY.find(" ") + 1;
				lineY = lineY.substr(0, pos);
				Xcoor[k - 3] = stod(lineX);
				Ycoor[k - 3] = stod(lineY);
				Zcoor[k - 3] = stod(lineZ);
			}
		}
		if (i >= fshots) //отсечение первых положений 
		{
			for (int j = 0; j < partnum; j++) //пробегание по всем парам частиц и внесение каждой пары в определённый слой
			{
				for (int l = j + 1; l < partnum; l++)
				{
					double r = 0; int num = 0;
					r = dist(Xcoor[j], Xcoor[l], Ycoor[j], Ycoor[l], Zcoor[j], Zcoor[l], length);
					num = floor(r / delta); 
					if (num < 100)
					{
						g[num] = g[num] + 1;
					}
				}
			}
		}
	}
	for (int i = 0; i < sloi; i++) //расчёт функции распределения и вывод в файл функции распределения
	{
		g[i] = g[i] / (shots - fshots);
		g[i] = raspred(g[i], i, delta, pairdensity);
		fout << ((i + 0.5) * delta) << " "<< g[i] << endl;
	}
	fout.close();
	fin.close();
	cout << pairdensity;
	system("pause");
	return (0);
}

