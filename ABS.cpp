#include <iostream>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <ctime>
using namespace std;


double* read(char x[10]) // функция прочтения файла конфигурации
{
	double conf[6]; 
	ifstream fin(x);
	fin >> conf[1] >> conf[2] >> conf[3] >> conf[4] >> conf [5];
	fin.close();
	return conf;
}


bool distcheck(double x1, double x2, double y1, double y2, double z1, double z2, double molrad, double length) //расчёт растояний между частицами и ответ на вопрос: не слишком ли близко прыгнула одна частица к другой
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

double jump(double coor,double dcoor, double length) //функция прыжка частицей, применяется покоординатно
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
	int n = ceil(pow(partnum, 1.0 / 3)); //создание сетки на кубе
	double griddist = length / (n);
	
	for (int i = 1; i <= partnum; i++) //расстановка частиц в начальное полежение (на сетку куба)
	{
		Xcoor[i] = griddist * ((i % n));
		Ycoor[i] = griddist * (((int) (i / n)) % n);
		Zcoor[i] = griddist * ((int) (i / (n*n)) % n);
	}
	ofstream fout;
	fout.open("1.XMOL");
	fout << partnum << endl; // вывод в файл начального положения
	fout << "STEP=" << 0 << endl;
	for (int i = 1; i <= partnum; i++)
	{
		fout << "Ar " << Xcoor[i] << " " << Ycoor[i] << " " << Zcoor[i] << endl;
	}
	
	for (int stepcount = 1; stepcount <= steps; stepcount++)
	{
		int r; 
		r = rand() % int (partnum) + 1; //определение частицы для прыжка
		

		double Xcoorrem, Ycoorrem, Zcoorrem;
		Xcoorrem = Xcoor[r]; Ycoorrem = Ycoor[r]; Zcoorrem = Zcoor[r]; //запоминание положения для возврата в случае отказа от прыжка
		
		Xcoor[r] = jump(Xcoor[r], dcoor, length);
		Ycoor[r] = jump(Ycoor[r], dcoor, length); //ïпрыжок частицей
		Zcoor[r] = jump(Zcoor[r], dcoor, length);
		
		bool t=true; //флаг принятия шага
		
		for (int k = 1; k <= partnum; k++) //проверка частиц на соприкосновение
		{	
		
			if (k == r)
			{
				t = true; 
			}
			else
			{
				t = distcheck(Xcoor[k], Xcoor[r], Ycoor[k], Ycoor[r], Zcoor[k], Zcoor[r], molrad, length); 
			}

			if (t == false)
			{ 
				break;
			}
		}

		if (t == false) //принятие или откат шага
		{
			badsteps = badsteps + 1;
			Xcoor[r] = Xcoorrem;
			Ycoor[r] = Ycoorrem;
			Zcoor[r] = Zcoorrem;
		}
		
		if ((stepcount % 125) == 0) //вывод в файл положений
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
