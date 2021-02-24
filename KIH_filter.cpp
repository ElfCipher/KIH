#include "pch.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <chrono>
//#include <xmmintrin.h>

const float PI = 3.141592653;
const uint32_t length = 100000;

class Timer
{
private:
	// Псевдонимы типов используются для удобного доступа к вложенным типам
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1>>;

	std::chrono::time_point<clock_t> m_beg;

public:
	Timer() : m_beg(clock_t::now())
	{
	}

	void reset()
	{
		m_beg = clock_t::now();
	}

	double elapsed() const
	{
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};

class KIH
{
public:

	float Fd; 			// Частота дискретизации, кГц
	float Fs; 			// Граничные частоты полосы пропускания, кГц
	float Fx; 			// Ширина полосы перехода, кГц
	uint16_t filter_order;

	KIH(float Fs, float Fx, float Fd)
	{
		this->Fs = Fs;
		this->Fx = Fx;
		this->Fd = Fd;
		this->filter_order = 53;
		h.resize(filter_order);
	}
	KIH (float Fs, float Fx, float Fd, uint16_t filter_order) : KIH(Fs, Fx, Fd)
	{
		this->filter_order = filter_order;
		h.resize(filter_order);
	}

	void GetFilteredSignal(std::vector <float> & const data_in, std::vector <float> &data_out)
	{
		// Применим метод взвешивания
		calculate_coefficients();

		for (size_t i = 0; i < length; i++)
		{
			float tmp = 0;

			for (size_t j = 0; j < filter_order && i >= j; j++)
				tmp += h[j] * data_in[i - j];

			data_out.push_back(tmp);
		}
	}

private:

	std::vector <float> h; 					// Импульсная характеристика фильтра

	void calculate_coefficients()
	{
		std::vector <float> h_D(filter_order, 0); 		// Идеальная импульсная характеристика
		std::vector <float> w(filter_order, 0);		// Весовая функция Хэмминга

		filter_order--;
		uint16_t n = filter_order / 2;
		float Fc = (Fs + Fx / 2) / Fd; // Центр полосы перехода
		// Для среднего коэффициента
		h_D[n] = 2 * Fc;
		w[n] = 0.54 + 0.46*cos(2 * PI * 0 / filter_order);
		h[n] = h_D[n] * w[n];

		// Для остальных
		for (size_t i = 1; i <= n; i++)
		{
			h_D[n-i] = sin(2 * PI*Fc*i) / (PI * i);
			h_D[n+i] = h_D[n-i]; 	// Так как значения коэффициентов симметричны относительно среднего элемента

			w[n-i] = 0.54 + 0.46*cos(2 * PI * i / filter_order);
			w[n+i] = w[n-i];

			h[n-i] = h_D[n-i] * w[n-i];
			h[n+i] = h[n-i];
		}
		filter_order++;
	}
	
};

void generate_signal(std::vector <float> & unfiltered_signal, float amplitude, float frequency, float Fd) // Записываем в один массив несколько синусов
{
	if (amplitude <= 0 || frequency <= 0)
	{
		printf("Must be positive! \n");
		return;
	}

	for (size_t i = 0; i < length; i++)
		unfiltered_signal[i] += amplitude * sin(2 * PI * frequency * i / Fd);
}

int main(void)
{
	float Fd = 8, 			// Частота дискретизации, кГц
		  Fs = 1.5, 		// Граничные частоты полосы пропускания, кГц
	      Fx = 0.5; 		// Ширина полосы перехода, кГц
	uint16_t filter_order = 23;

	std::vector <float> unfiltered_signal(length, 0);
	std::vector <float> filtered_signal;
	// Три сигнала с разными частотой и амплитудой
	generate_signal(unfiltered_signal, 10, 1, Fd); // будет пропускать этот сигнал
	generate_signal(unfiltered_signal, 2, 2, Fd);
	generate_signal(unfiltered_signal, 1, 3, Fd);

	Timer t;

	KIH LPF(Fs, Fx, Fd, filter_order);
	LPF.GetFilteredSignal(unfiltered_signal, filtered_signal);

	std::cout << "Time elapsed: " << t.elapsed() << '\n';

	return 0;
}
