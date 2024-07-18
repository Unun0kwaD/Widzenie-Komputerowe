# Compilation

arm-linux-gnueabihf-gcc -mfpu=neon -static -ftree-vectorize main.c -o main -lm

# Example usage

qemu-arm ./main d 10 west_2.ppm west_dilation2.pgm

qemu-arm ./zad1 e 10 west_2.ppm west_erosion2.pgm
