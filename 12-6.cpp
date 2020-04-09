#include <iostream>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <ctime>
using namespace std;

double* read(char x[10]) //функция для прочтения конфигурационного файла
{
	double conf[5];
	ifstream fin(x);
	fin >> conf[0] >> conf[1] >> conf[2] >> conf[3] >> conf[4];
	fin.close();
	return conf;
}

double dist(double x1, double x2, double y1, double y2, double z1, double z2, double length) //функция расчёта дистанции между частицами методом наиближайшего образа
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
	return (r);
}

double jump(double coor, double dcoor, double length) // функция прыжка для выбранной частицы (для покоординатное использование)
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

double U(double r, double epsilon_T) //расчёт энергии между 2-мя частицами
{
	double t = 1 / r;
	double U = 4 * epsilon_T* (pow(t, 12) - pow(t, 6));
	return (U);
}

bool jumpcheck(double dE) //проверка прыжка на подтверждение (Если энергия новой конфигурации меньше, то прыжок гарантированно принят, иначе считаем вероятность и кидаем случайное число для принития или непринятия)
{
	double pjump; bool t;
	pjump = exp(-dE);
	if (pjump > 1)
	{
		pjump = 1;
	}
	double p;
	p = rand(); p = p / RAND_MAX;
	
	if (p < pjump)
	{
		t = true;
	}
	else
	{
		t = false;
	}

	return (t);
}

int main()
{
	srand(static_cast<unsigned int>(time(0)));
	char way[10]; double partnum, steps, length, dcoor, molrad; int csteps = 0; int badsteps = 0; int shots = 0; double Eavrshot = 0;
	double Esum = 0; double sigma = 3.4 * pow(10, -10); double k = 1.38 * pow(10, -23); double T_norm = 120; double epsilon = T_norm * k; double Eusr = 0; double T = 0;
	cout << "cway = ";
	cin >> way;
	cout << endl << "T = ";
	cin >> T;
	double epsilon_T = T_norm / T;
	partnum = read(way)[0];
	steps = read(way)[1];
	length = read(way)[2];
	dcoor = read(way)[3];
	molrad = read(way)[4];
	double* Xcoor = new double[partnum];
	double* Ycoor = new double[partnum];
	double* Zcoor = new double[partnum];
	
	int n = ceil(pow(partnum, 1.0 / 3)); // количество узлов в сетке на грань куба
	
	double griddist = length / (n);

	for (int i = 0; i < partnum; i++) //расстановка по начальной сетке
	{
		Xcoor[i] = griddist * ((i % n));
		Ycoor[i] = griddist * (((int)(i / n)) % n);
		Zcoor[i] = griddist * ((int)(i / (n * n)) % n);
	}

	ofstream Eout; //файл энергий
	Eout.open("E-1.txt");
	ofstream fout; //файл положений
	fout.open("6-12-1.XMOL");
	fout << partnum << endl;
	fout << "STEP=" << 0 << endl;
	for (int i = 0; i < partnum; i++)
	{
		fout << "Ar " << Xcoor[i] << " " << Ycoor[i] << " " << Zcoor[i] << endl;
	}

	for (int i=0; i<partnum; i++) //расчёт суммарной энергии начальной конфигурации
	{
		for (int j = i + 1; j < partnum; j++)
		{
			double d = 0;
			d = dist(Xcoor[i], Xcoor[j], Ycoor[i], Ycoor[j], Zcoor[i], Zcoor[j], length);
			Esum = Esum + U(d,epsilon_T);
		}
	}

	Eout << 0 << " " << Esum << endl;

	for (int stepcount = 0; stepcount < steps; stepcount++)
	{
		int l;
		l = rand() % int(partnum); //выбрали частицу

		double Xcoorrem, Ycoorrem, Zcoorrem;
		Xcoorrem = Xcoor[l]; Ycoorrem = Ycoor[l]; Zcoorrem = Zcoor[l]; //запомнили положение частицы

		Xcoor[l] = jump(Xcoor[l], dcoor, length);
		Ycoor[l] = jump(Ycoor[l], dcoor, length); //прыжок частицей
		Zcoor[l] = jump(Zcoor[l], dcoor, length);

		bool t = true; //флаг принятия шага
		
		double dEsum = 0; //измениение энергии от положения к положению
		for (int k = 0; k < partnum; k++)
		{
			if (k != l)
			{
				double d_before = dist(Xcoorrem, Xcoor[k], Ycoorrem, Ycoor[k], Zcoorrem, Zcoor[k], length);
				double d_after = dist(Xcoor[l], Xcoor[k], Ycoor[l], Ycoor[k], Zcoor[l], Zcoor[k], length);
				double U1 = U(d_before, epsilon_T);
				double U2 = U(d_after, epsilon_T);
				dEsum = dEsum + (U2-U1); 
			}
		}

		t = jumpcheck(dEsum);

		if (t == false) //принятие или непринятие шага
		{
			Xcoor[l] = Xcoorrem;
			Ycoor[l] = Ycoorrem;
			Zcoor[l] = Zcoorrem;
			badsteps = badsteps +1;
		}
		else
		{
			Esum = Esum + dEsum;
		}
		


		if 
			(((stepcount+1) % 1000) == 1) //вывод в файл положений раз в 1000 шагов
		{
			fout << partnum << endl;
			fout << "STEP=" << stepcount+1 << endl;
			shots = shots + 1;
			for (int i = 0; i < partnum; i++)
			{
				fout << "Ar " << Xcoor[i] << " " << Ycoor[i] << " " << Zcoor[i] << endl;
			}
		} 

		if (stepcount % 10000 == 0) //вывод шага для удобства пользования
		{
			cout << stepcount << endl;
		}

		if (stepcount % 2500 == 0) // перепроверка суммы
		{
			double Esumcheck = 0; double otn = 0;
			for (int i = 0; i < partnum; i++)
			{
				for (int j = i + 1; j < partnum; j++)
				{
					Esumcheck = Esumcheck + U(dist(Xcoor[i], Xcoor[j], Ycoor[i], Ycoor[j], Zcoor[i], Zcoor[j], length), epsilon_T);
				}
			}


			otn = Esumcheck / Esum;
			
			
			if ((otn < 0.9999) || (otn > 1.0001))
			{
				cout << "!!!bad sum!!!   step = " << stepcount << "   otn = " << otn << "   Esumcheck = " << Esumcheck << "   Esum = " << Esum << "    dE="<< dEsum  << endl;

				Esum = Esumcheck;
			}
		}

		if (stepcount > 50000) // средняя энергия
		{
			Eusr = Eusr + Esum;
			csteps = csteps + 1;
			if (stepcount % 500 == 0)
			{
				Eout << stepcount + 1 << "E= " << Esum << " dE=" << dEsum << endl;// вывод энергий в файл
			}
		}
	}
	Eout.close();
	fout.close();
	fout.open("c1.txt");
	fout << partnum << " " << steps << " " << length << " " << dcoor << " " << molrad << " " << shots << endl;
	fout.close();
	cout << endl << partnum << " " << steps << " " << length << " " << dcoor << " " << molrad << endl;
	cout << "E=" << Eusr/csteps << endl;
	cout << badsteps / steps;
	system("pause");
	return 0;
}
