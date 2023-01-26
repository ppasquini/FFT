compression: main.cpp FFT_2D.cpp FFT_2D.hpp IMAGE_COMPRESSION.cpp IMAGE_COMPRESSION.hpp FFT_traits.hpp
	g++ main.cpp FFT_2D.cpp IMAGE_COMPRESSION.cpp -fopenmp -o compression -I ${mkEigenInc}
	mkdir -p Matrix

FFT: main.cpp FFT_2D.cpp FFT_2D.hpp FFT_traits.hpp
	g++ main_FFT_2D.cpp FFT_2D.cpp -fopenmp -o FFT -I ${mkEigenInc}