/* $Id$ */
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#define MIN(a,b)	((a) < (b) ? (a) : (b))
#define MAX(a,b)	((a) < (b) ? (b) : (a))

typedef struct __attribute__((__packed__)) {
	uint8_t		bfType1;
	uint8_t		bfType2;
	uint32_t	bfSize;       
	uint16_t	bfReserved1;      
	uint16_t	bfReserved2;          
	uint32_t	bfOffBits;
} FILEHEADER;

typedef struct __attribute__((__packed__)) {
	uint32_t	biSize;      
	uint32_t	biWidth;
	uint32_t	biHeight;
	uint16_t	biPlanes;        
	uint16_t	biBitCount;        
	uint32_t	biCompression;        
	uint32_t	biSizeImage;       
	uint32_t	biXPelsPerMeter; 
	uint32_t	biYPelsPerMeter;       
	uint32_t	biClrUsed;        
	uint32_t	biClrImportant;    
} INFOHEADER;

typedef struct __attribute__((__packed__)) {
	uint8_t	b;
	uint8_t	g;
	uint8_t	r;
} IMAGE;

typedef struct {
	FILEHEADER fh;
	INFOHEADER ih;
	IMAGE *pix;
	FILE *dest;
	int o_w;
	int o_h;
	int n_w;
	int n_h;
	int buf;
} DATA;

IMAGE *getData(FILE *, FILEHEADER, INFOHEADER);
IMAGE *makeImage(int, int, int);
FILEHEADER getFileHeader(FILE *);
INFOHEADER getInfoHeader(FILE *);
void printFileHeader(FILEHEADER);
void printInfoHeader(INFOHEADER);
void printHex(IMAGE *, int, int);
void resize(char *, char *, int, int, int, void (*)(FILEHEADER, INFOHEADER, IMAGE *, FILE *, int, int, int, int, int));

int main(int argc, char *argv[]) {
	void nearestNeighbour(FILEHEADER, INFOHEADER, IMAGE *, FILE *, int, int, int, int, int);	
	void bilinearInterpolation(FILEHEADER, INFOHEADER, IMAGE *, FILE *, int, int, int, int, int);
	void bicubicInterpolation(FILEHEADER, INFOHEADER, IMAGE *, FILE *, int, int, int, int, int);
	
	int c, buf = 1, time = 0, w = 0, h = 0, s = 0;
	clock_t start, end;
	char *src = NULL, *dest = NULL;
	void (*func)(FILEHEADER, INFOHEADER, IMAGE *, FILE *, int, int, int, int, int) = bilinearInterpolation;	
	//src = argv[1];
	//dest = argv[2];
	while (optind < argc) {
		if ((c = getopt(argc, argv, "i:d:btw:h:s:c")) != -1) {
			switch (c) {
				case 'i': { 
					FILE *f;
					if ((f = fopen(optarg, "rb")) == NULL) {
						printf("Bad source file: %s\n", optarg);
						return -1;
					}
					FILEHEADER fh = getFileHeader(f);
					INFOHEADER ih = getInfoHeader(f);
					printFileHeader(fh);
					printInfoHeader(ih);
					return 0;
				}
				case 'd': {
					FILE *f;
					if ((f = fopen(optarg, "rb")) == NULL) {
						printf("Bad source file: %s\n", optarg);
						return -1;
					}
					FILEHEADER fh = getFileHeader(f);
					INFOHEADER ih = getInfoHeader(f);
					IMAGE *img = getData(f, fh, ih);
					printHex(img, ih.biWidth, ih.biHeight);
					return 0;
				}
				case 'b':
					buf = 0; // Turns OFF buffering
					break;
				case 't':
					time = 1;
					break;
				case 'w':
					w = (int)atol(optarg);
					break;
				case 'h':
					h = (int)atol(optarg);
					break;
				case 's':
					s = (int)atol(optarg);
					break;
				case 'c':
					printf("Select scaling algorithm:\n(1) Nearest Neighbour\n(2) Bilinear Interpolation\n(3) Bicubic Interpolation\n");
					int input = 0;
					if (scanf("%d", &input) != 1)
						printf("Invalid input\n");
					switch (input) {
						case 1:
							func = nearestNeighbour;
							break;
						case 2:
							func = bilinearInterpolation;
							break;
						case 3:
							func = bicubicInterpolation;
						default:
							break;
					}
				default:
					break;
			}
		}
		else {
			if (src == NULL)
				src = argv[optind];
			else if (dest == NULL)
				dest = argv[optind];
			else {
				printf("Usage: %s\n", argv[optind]);
				return -1;
			}
			optind++;
		}
	}

	if (s) 
		w = h = s*-1;
	else if (h == 0 || w == 0) {
		fprintf(stderr, "ERROR: Invalid width, height or scale values.\n");
		return -1;
	}
	start = clock();
	resize(src, dest, w, h, buf, func);
	end = clock();
	if (time)
		printf("Took %.3fms\n", (double)(end-start)/(CLOCKS_PER_SEC/1000));
	return 0;
}

void resize(char *s_src, char *s_dest, int n_w, int n_h, int buf, void (*func)(FILEHEADER, INFOHEADER, IMAGE *, FILE *, int, int, int, int, int)) {	
	FILE *src, *dest;
	if ((src = fopen(s_src, "rb")) == NULL) {
		printf("Bad source file: %s\n", s_src);
		return;
	}
	if ((dest = fopen(s_dest, "wb")) == NULL)  {
		printf("Bad destination file: %s\n", s_dest);
		return;
	}
	FILEHEADER fh = getFileHeader(src);
	INFOHEADER ih = getInfoHeader(src);
	IMAGE *pix = getData(src, fh, ih);
	int o_w = ih.biWidth, o_h = ih.biHeight;
	if (n_w > 0 && n_h > 0) {
		ih.biWidth = n_w;
		ih.biHeight = n_h;
	}
	else if (n_w != n_h) {
		printf("Invalid scale or dimensions\n");
		return;
	}
	else {
		ih.biWidth *= ((double)n_w/100)*-1;
		ih.biHeight *= ((double)n_w/100)*-1;
		n_w = ih.biWidth;
		n_h = ih.biHeight;
	}
	fwrite(&fh, 1, sizeof(FILEHEADER), dest);
	fwrite(&ih, 1, sizeof(INFOHEADER), dest);

	func(fh, ih, pix, dest, o_w, o_h, n_w, n_h, buf);

	fclose(src);
	fclose(dest);
}

void nearestNeighbour(FILEHEADER fh, INFOHEADER ih, IMAGE *pix, FILE *dest, int o_w, int o_h, int n_w, int n_h, int buf) {	
	double xrat = (double)(o_w-1)/(double)n_w;
	double yrat = (double)(o_h-1)/(double)n_h;
	int x, y, i, j;
	switch (ih.biBitCount) {
		case 16: {
			int pad = (2*n_w) & 3;
			if (pad)
				pad = 4-pad;
			uint16_t *buffer;
			if (buf)
				buffer = malloc(2*n_w);
			for (i = 0; i <n_h; i++) {
				for (j = 0; j < n_w; j++) {
					x = (int)(j*xrat);
					y = (int)(i*yrat);
					uint16_t container = 0;
					container |= ((pix+(y*o_w)+x)->b) << 11;
					container |= ((pix+(y*o_w)+x)->g) << 5;
					container |= ((pix+(y*o_w)+x)->r);
					if (buf)
						*(buffer+j) = container;
					else
						fwrite(&container, 1, 2, dest);
				}
				if (buf)
					fwrite(buffer, 1, 2*n_w, dest);
				fwrite(&pad, 1, pad, dest); 
			}
			break;
		}
		case 24: { 
			int pad = (sizeof(IMAGE)*n_w) & 3;
			if (pad)
				pad = 4-pad;
			IMAGE *buffer;
			if (buf)
				buffer = malloc(sizeof(IMAGE)*n_w);
			for (i = 0; i < n_h; i++) {
				for (j = 0; j < n_w; j++) {
					x = (int)(j*xrat);
					y = (int)(i*yrat);
					if (buf)
						*(buffer+j) = *(pix+(y*o_w)+x);
					else
						fwrite(pix+(y*o_w)+x, 1, sizeof(IMAGE), dest);
				}
				if (buf)
					fwrite(buffer, 1, sizeof(IMAGE)*n_w, dest);
				fwrite(&pad, 1, pad, dest);	
			}
			break;
		}
		default:
			printf("Unsupported: %dbit\n", ih.biBitCount);
			return;
	}
}

void bilinearInterpolation(FILEHEADER fh, INFOHEADER ih, IMAGE *pix, FILE *dest, int o_w, int o_h, int n_w, int n_h, int buf) {
	#define getelm(x)	((((pix+index)->x)*(1-xdiff)*(1-ydiff))+(((pix+index+1)->x)*(xdiff)*(1-ydiff))+(((pix+index+o_w)->x)*(1-xdiff)*(ydiff))+((pix+index+o_w+1)->x)*xdiff*ydiff)
	double xrat = (double)(o_w-1)/(double)n_w;
	double yrat = (double)(o_h-1)/(double)n_h;
	int x, y, i, j, index;
	double xdiff, ydiff;
	switch(ih.biBitCount) {
		case 16: {
			int pad = (2*n_w) & 3;
			if (pad)
				pad = 4-pad;
			uint16_t *buffer;
			if (buf)
				buffer = malloc(2*n_w);
			for (i = 0; i < n_h; i++) {
				for (j = 0; j < n_w; j++) {
					x = (int)(j*xrat);
					y = (int)(i*yrat);
					xdiff = (xrat*j)-x;
					ydiff = (yrat*i)-y;
					index = y*o_w+x;
					uint16_t container = 0;
					double bc = getelm(b);
					double gc = getelm(g);
					double rc = getelm(r);
					if (bc > 31) {
						printf("b overflow\n");
						bc = 15;
					}
					if (gc > 63) {
						printf("g overflow\n");
						gc = 15;
					}
					if (rc > 31) {
						printf("r overflow\n");
						rc = 15;
					}
					container |= (int)bc << 11;
					container |= (int)gc << 5;
					container |= (int)rc;
					if (buf)
						*(buffer+j) = container;
					else
						fwrite(&container, 1, 2, dest);
				}
				if (buf)
					fwrite(buffer, 1, 2*n_w, dest);
				fwrite(&pad, 1, pad, dest);
			}
			break;
		}
		case 24: {
			int pad = (sizeof(IMAGE)*n_w) & 3;
			if (pad)
				pad = 4-pad;
			IMAGE *buffer;
			if (buf)
				buffer = malloc(sizeof(IMAGE)*n_w);
			for (i = 0; i < n_h; i++) {
				for (j = 0; j < n_w; j++) {
					x = (int) (xrat*j);
					y = (int) (yrat*i);
					xdiff = (xrat*j)-x;
					ydiff = (yrat*i)-y;
					index = y*o_w+x;
					(buffer+j)->r = getelm(r);
					(buffer+j)->g = getelm(g);
					(buffer+j)->b = getelm(b);
				}
				if (buf)
					fwrite(buffer, 1, sizeof(IMAGE)*n_w, dest);
				fwrite(&pad, 1, pad, dest);
			}
			break;
		}
		default:
			break;
	}
}

void bicubicInterpolation(FILEHEADER fh, INFOHEADER ih, IMAGE *pix, FILE *dest, int o_w, int o_h, int n_w, int n_h, int buf) {
#define posPix(X)	((y-1+jj)*o_w + (x+(X)) + k)
	
	double xrat = (double)(o_w)/(double)n_w;
	double yrat = (double)(o_h)/(double)n_h;
	int x, y, i, j, k, jj, z, index;
 	int a1[3], a2[3], a3[3], C[5][3] = {0,0,0};
	double xdiff, ydiff;

	int pad = (sizeof(IMAGE)*n_w) & 3;
	if (pad)
		pad = 4-pad;
	IMAGE *buffer = malloc(sizeof(IMAGE)*n_w*n_h); // Bicubic Interpolation can't do line buffers, it handles chunks
	for (i = 0; i < n_h; i++) {
		for (j = 0; j < n_w; j++) {
			x = (int) (xrat*j);
			y = (int) (yrat*i);
			xdiff = (xrat*j)-x;
			ydiff = (yrat*i)-y;
			index = y*o_w+x;

			for (k = 0; k < 3; k++) {
				for (jj = 0; jj < 4; jj++) {
					int a0[3] = {((pix+posPix(0))->r), ((pix+posPix(0))->g), ((pix+posPix(0))->b)}; 

					int d0[3] = {((pix+posPix(-1))->r - a0[0]), ((pix+posPix(-1))->g - a0[1]), ((pix+posPix(-1))->b - a0[2])};
					int d2[3] = {((pix+posPix(1))->r - a0[0]), ((pix+posPix(1))->g - a0[1]), ((pix+posPix(1))->b - a0[2])};
					int d3[3] = {((pix+posPix(2))->r - a0[0]), ((pix+posPix(2))->g - a0[1]), ((pix+posPix(2))->b - a0[2])};

					for (z=0; z<3; z++)
						a1[z] = (-1.0/3*d0[z]+d2[z] -1.0/6*d3[z]);
					for (z=0; z<3; z++)
						a2[z] = (1.0/2*d0[z]+1.0/2*d2[z]);
					for (z=0; z<3; z++)
						a3[z] = (-1.0/6*d0[z]-1.0/2*d2[z]+1.0/6*d3[z]);

					for (z=0; z<3; z++)
						C[jj][z] = (a0[z]+a1[z]*xdiff+a2[z]*xdiff*xdiff+a3[z]*xdiff*xdiff*xdiff);

					for (z=0; z<3; z++)
						d0[z] = C[0][z] - C[1][z];
					for (z=0; z<3; z++)
						d0[z] = C[2][z] - C[1][z];
					for (z=0; z<3; z++)
						d0[z] = C[3][z] - C[1][z];				
				
					for (z=0; z<3; z++)
						a0[z] = C[1][z];
				
					for (z=0; z<3; z++)
						a1[z] = -1.0/3*d0[z] +d2[z] -1.0/6*d3[z];
					for (z=0; z<3; z++)
						a2[z] = 1.0/2*d0[z] + 1.0/2*d2[z];
					for (z=0; z<3; z++)
						a3[z] = -1.0/6*d0[z] - 1.0/2*d2[z] +1.0/6*d3[z];

					(buffer+(i*n_w+j+k))->r = (a0[0]+a1[0]*ydiff+a2[0]*ydiff*ydiff+a3[0]*ydiff*ydiff*ydiff);
					(buffer+(i*n_w+j+k))->g = (a0[1]+a1[1]*ydiff+a2[1]*ydiff*ydiff+a3[1]*ydiff*ydiff*ydiff);
					(buffer+(i*n_w+j+k))->b = (a0[2]+a1[2]*ydiff+a2[2]*ydiff*ydiff+a3[2]*ydiff*ydiff*ydiff);
				}
			}
		}
	}
	for (z=0; z<n_h; z++) {
		fwrite((buffer+(z*n_w)), 1, sizeof(IMAGE)*n_w, dest);
		fwrite(&pad, 1, pad, dest);
	}
}


FILEHEADER getFileHeader(FILE *f) {
	FILEHEADER fh;
	long loc = ftell(f);
	fseek(f, 0, 0);
	fread(&fh, 1, sizeof(FILEHEADER), f);
	fseek(f, loc, 0);
	return fh;
}

INFOHEADER getInfoHeader(FILE *f) {
	INFOHEADER ih;
	long loc = ftell(f);
	fseek(f, sizeof(FILEHEADER), 0);
	fread(&ih, 1, sizeof(INFOHEADER), f);
	fseek(f, loc, 0);
	return ih;
}

IMAGE *getData(FILE *f, FILEHEADER fh, INFOHEADER ih) {
	long loc = ftell(f);
	fseek(f, fh.bfOffBits, 0);
	IMAGE *im = malloc(sizeof(IMAGE)*ih.biWidth*ih.biHeight);
	int i, j;
	switch (ih.biBitCount) {
		case 16: {
			uint16_t container;
			for (i = 0; i < ih.biHeight;  i++)
				for (j = 0; j < ih.biWidth; j++) {
					fread(&container, 1, sizeof(container), f);
					(im+((i*ih.biWidth)+j))->r = container & 31;
					(im+((i*ih.biWidth)+j))->g = (container & 2016) >> 5;
					(im+((i*ih.biWidth)+j))->b = (container & 63488) >> 11;
				}
			break;
		}
		case 24:
			for (i = 0; i < ih.biHeight;  i++)
				for (j = 0; j < ih.biWidth; j++)
					fread((im+((i*ih.biWidth)+j)), 1, sizeof(IMAGE), f);
			break;
		default:
			printf("Unsupported file format: %dbit\n", ih.biBitCount);
			break;
	}
	fseek(f, loc, 0);
	return im;
}

IMAGE *makeImage(int r, int g, int b) {
	IMAGE *img = malloc(sizeof(IMAGE));
	img->r = r;
	img->g = g;
	img->b = b;
	return img;
}

void printHex(IMAGE *im, int width, int height) {
	int i, p = width*height;
	for (i = 0; i < p; i++)
			printf("%.2x%.2x%.2x%s", (im+i)->r, (im+i)->g, (im+i)->b, (((i+1)%width)==0)?"\n":" ");
}

void printFileHeader(FILEHEADER fh) {
	printf("\n-= Info Header =-\nbfType: %c%c, bfSize: %d, bfReserved1: %d, bfReserved2: %d, bfoffBits: %d.\n", \
		fh.bfType1, fh.bfType2, fh.bfSize, fh.bfReserved1, fh.bfReserved2, fh.bfOffBits);
}

void printInfoHeader(INFOHEADER ih) {
	printf("\n-= Image Header =-\nbiSize: %d, biWidth: %d, biHeight: %d, biPlanes: %d, biBitCount: %d, biBitCompr"\
		"ession: %d, biSizeImage: %d, biXPelsPerMeter: %d, biYPelsPerMeter: %d, biClrUsed: %d, biClrImportant: %d\n\n",\
		ih.biSize, ih.biWidth, ih.biHeight, ih.biPlanes, ih.biBitCount, ih.biCompression, ih.biSizeImage, ih.biXPelsPerMeter,\
		ih.biYPelsPerMeter, ih.biClrUsed, ih.biClrImportant);
}
