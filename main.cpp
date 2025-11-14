/*
 * main.cpp
 *
 *  Created on: 12 нояб. 2025 г.
 *      Author: Andrew
 */


#include "math.h"

#include "string.h"

#include <string>

#include <stdio.h>
#include <stdlib.h>

#include  <array>
#include  <vector>

extern "C"
{
#include "libbmp.h"
#include "cJSON.h"
#include "font.h"
#include "csv.h"
}


extern const unsigned char font8x8_ib8x8u_full[][8];



typedef struct
{
	double x1;
	double y1;
	double x2;
	double y2;
}
Vector_t;


class Vessel
{
	unsigned int id;
	double heading_degrees;
	double speed;
	double coords_x;
	double coords_y;

	double CPA;
	double TCPA;

public:

	Vessel(	unsigned int _id,
			double _coords_x,
			double _coords_y,
			double _speed,
			double _heading_degrees
		   )
	{
		id = _id;
		heading_degrees = _heading_degrees;
		speed = _speed;
		coords_x = _coords_x;
		coords_y = _coords_y;
		CPA = 0;
		TCPA = 0;
	}


	Vessel()
	{
		id = 0;
		heading_degrees = 0;
		speed = 0;
		coords_x = 0;
		coords_y = 0;
		CPA = 0;
		TCPA = 0;
	}

	int getID(){return id;}
	double getX(){return coords_x;}
	double getY(){return coords_y;}
	double getCourse(){return heading_degrees;}
	double getSpeed(){return speed;}

	void setCPA(double cpa)
	{
		CPA = cpa;
	}

	void setTCPA(double tcpa)
	{
		TCPA = tcpa;
	}


	double getCPA(void)
	{
		return CPA;
	}

	double getTCPA(void)
	{
		return TCPA;
	}

};


class Vector
{
private:

	//Decart system
	double x1;
	double y1;

	double x2;
	double y2;

	double length;
	double theta;

public:

	Vector(double in_x1, double in_y1, double in_x2, double in_y2)
	{
		x1 = in_x1;
		x2 = in_x2;
		y1 = in_y1;
		y2 = in_y2;

		double dx = x2-x1;
		double dy = y2-y1;

		length = sqrt(dx * dx + dy * dy);
		theta =  atan2(dy, dx);
	}

	double getLength()
	{
		return length;
	}


	static Vector VectorHL(double in_x, double in_y, double angle, double length)
	{
		double x1 = in_x;
		double y1 = in_y;

		double x2 = in_x + length * cos(angle);
		double y2 = in_y + length * sin(angle);

		return Vector(x1,y1,x2,y2);
	}



	Vector operator+(Vector &v2)
	{
		return Vector(x1+v2.x1, y1+v2.y1, x2+v2.x2, y2+v2.y2);
	};

	Vector operator-(Vector &v2)
	{
		return Vector(x1-v2.x1, y1-v2.y1, x2-v2.x2, y2-v2.y2);
	};

	double operator*(Vector &v2) //vector-vector scalar dot product
	{
		double lx1 = (x2 - x1);
		double ly1 = (y2 - y1);

		double lx2 = (v2.x2 - v2.x1);
		double ly2 = (v2.y2 - v2.y1);

		double l1 = sqrt(lx1*lx1 + ly1*ly1);
		double l2 = sqrt(lx2*lx2 + ly2*ly2);

		double angle1 = atan2(ly1, lx1);
		double angle2 = atan2(ly2, lx2);

		double theta = angle2 - angle1;

		return l1*l2*cos(theta);
	};

	Vector operator*(double k)  //vector-scalar dot product
	{
		double lx = k * (x2 - x1);
		double ly = k * (y2 - y1);

		return Vector(x1, y1, (x1 + lx), (y1 + ly));
	};

};





Vessel ownVessel;
//Vessel_t* targetVessels;

int targetVesselsBufferSize = 10;
int targetVesselsCount;

class BMP_Image
{
	bmp_img img;


public:

	void inline __attribute__((always_inline)) draw_pixel(int row, int col, unsigned int rgb);
	void draw_line(int row1, int col1, int row2, int col2, unsigned int rgb);
	void draw_circle(int row, int col, float radius, unsigned int rgb);
	void draw_symbol(uint8_t* drawTable, char sym, int row, int col, unsigned int rgb);
	void draw_text(const char* txt, int row, int col, unsigned int rgb);
	void img_write(const char* filename);

	BMP_Image(char* filename)
	{
		bmp_img_read(&img, filename);
	}

	BMP_Image(int width, int height)
	{
		bmp_img_init_df(&img, width, height);

		for(int row=0; row<img.img_header.biHeight; ++row)
		{
			for(int col=0; col<img.img_header.biWidth; ++col)
			{
				draw_pixel(row, col, 0x00);
			}
		}
	}

	int getHeight(void){ return img.img_header.biHeight;}
	int getWidth(void){ return img.img_header.biWidth;}
};



std::vector<Vessel> targetVessels;

int main(int argc, char** argv)
{
	if(argc < 2)
	{
		printf("No input file.");
		return 1;
	}



	cJSON* json_report;

	json_report = cJSON_CreateObject();

	FILE* file = fopen(argv[1], "r");

	int ownVesselFound = 0;

	char line[1024];
	if (file != NULL)
	{
		char** vessels;
		while (fgets(line, sizeof(line), file))
		{
			 vessels = parse_csv(line);

			 int c = 0;
			 char* p = vessels[0];
			 while(p)
			 {
				 c++;
				 p = vessels[c];
			 }

			 if(c < 5)
			 {
				 printf("Error parsing CSV.\r\n");
				 return 1;
			 }



			 Vessel vessel = Vessel (  std::stoi( std::string(vessels[0]) ),
									   std::stod( std::string(vessels[1]) ),
									   std::stod( std::string(vessels[2]) ),
									   std::stod( std::string(vessels[3]) ),
									   std::stod( std::string(vessels[4]) ) );


			 if(vessel.getID() == 0)
			 {
				 ownVessel = vessel;

				 ownVesselFound = 1;
			 }
			 else
			 {
				 targetVessels.push_back(vessel);
			 }

			 free_csv_line(vessels);
			 vessels = NULL;
		}


		fclose(file);
	}
	else
	{
		printf("Cannot open file.\r\n");
		return 1;
	}

	if(!ownVesselFound)
	{
		printf("Own vessel not found.\r\n");
		return 1;
	}

	//create output bmp
	BMP_Image bmp_output = BMP_Image(600,600);


	//Draw own vessel
	bmp_output.draw_circle(bmp_output.getHeight() / 2, bmp_output.getWidth() / 2, 4, 0xFFFFFF);

	//Draw own vessel course line
	bmp_output.draw_line(bmp_output.getHeight() / 2, bmp_output.getWidth() / 2,
				0,
				bmp_output.getHeight() / 2,
				0x555555);

	//Draw own vessel speed line
	bmp_output.draw_line(bmp_output.getHeight() / 2,  bmp_output.getWidth() / 2,
				bmp_output.getHeight()  / 2 - ownVessel.getSpeed(),
				bmp_output.getHeight()  / 2,
				0xFFFFFF);

	Vector vector_speed1 = Vector::VectorHL ( 0,
											 0,
											 (90.0 - ownVessel.getCourse()) * (M_PI / 180.0F),
											 ownVessel.getSpeed() );

	for(unsigned int v = 0; v<targetVessels.size(); ++v)
	{
		Vector vector_speed2 = Vector::VectorHL (0,
												0,
												(90.0F - targetVessels.at(v).getCourse()) * (M_PI / 180.0F),
												targetVessels.at(v).getSpeed() );

		Vector vector_v = vector_speed2 - vector_speed1;
		//Vector vector_r = Vector(targetVessels.at(v) .getX(), targetVessels.at(v).getY(), ownVessel .getX(), ownVessel.getY());
		Vector vector_r = Vector(ownVessel.getX(), ownVessel.getY(), targetVessels.at(v).getX(), targetVessels.at(v).getY());


		double tcpa = -(vector_v * vector_r) / (vector_v.getLength() * vector_v.getLength());

		Vector t = (vector_v * tcpa);

		double cpa = (vector_r + t).getLength();

		targetVessels.at(v).setCPA(cpa);
		targetVessels.at(v).setTCPA(tcpa);



		//Check for collision risk.
		uint32_t vesselColor = 0x0000FF;
		uint32_t courseColor = 0x000066;
		uint32_t cpaColor = 0x00FFFF;
		if((tcpa > 0.0) && (tcpa < 30.0) && (cpa < 50.0))
		{
			vesselColor = 0xFF0000;
			courseColor = 0x660000;
			cpaColor = 0xFFFF00;
		}

		int cpaColOwn;
		int cpaRowOwn;

		int cpaColTarget;
		int cpaRowTarget;

		//Calculate own vessel CPA point
		Vector vCPA = (vector_speed1 * tcpa);

		cpaColOwn =  (int)ownVessel.getX() + bmp_output.getWidth()/2;
		cpaRowOwn =  bmp_output.getHeight() / 2 - (vCPA.getLength() + (int)ownVessel.getY());

		//Draw own vessel CPA point
		if(tcpa >= 0.0)
		{
			bmp_output.draw_circle(cpaRowOwn, cpaColOwn, 2, cpaColor);
		}

		//Calculate target vessel CPA point
		vCPA = (vector_speed2 * tcpa);
		cpaColTarget =  (int)targetVessels.at(v) .getX() + vCPA.getLength()*cos( (90.0 - targetVessels.at(v).getCourse()) * (M_PI/180.0)) +  bmp_output.getWidth()/2;
		cpaRowTarget =  bmp_output.getHeight()/2 - ((int)targetVessels.at(v).getY() + vCPA.getLength()*sin( (90.0F - targetVessels.at(v).getCourse()) * (M_PI/180.0)));

		//Draw target vessel CPA point
		if(tcpa >= 0.0)
		{
			bmp_output.draw_circle(cpaRowTarget, cpaColTarget, 2, cpaColor);
		}

		//Draw target vessel
		int vesselRow =  bmp_output.getHeight() / 2 - targetVessels.at(v).getY();
		int vesselCol =  bmp_output.getWidth() / 2  + targetVessels.at(v).getX();
		bmp_output.draw_circle(vesselRow, vesselCol, 4, vesselColor);

		//Set course line len
		int lineLen = sqrt( bmp_output.getHeight()* bmp_output.getHeight() +
					   bmp_output.getWidth()* bmp_output.getWidth());

		//Draw target vessel course line
		bmp_output.draw_line(vesselRow, vesselCol,
				vesselRow - (int)(lineLen * sin( (90.0F - targetVessels.at(v).getCourse()) * (M_PI/180.0)) ),
				vesselCol + (int)(lineLen * cos( (90.0F - targetVessels.at(v).getCourse()) * (M_PI/180.0)) ),
				courseColor);

		//Draw target vessel speed line
		bmp_output.draw_line(vesselRow, vesselCol,
				vesselRow - (int)(targetVessels.at(v).getSpeed() * sin( (90.0F - targetVessels.at(v).getCourse()) * (M_PI/180.0)) ),
				vesselCol + (int)(targetVessels.at(v).getSpeed() * cos( (90.0F - targetVessels.at(v).getCourse()) * (M_PI/180.0)) ),
				vesselColor);


		if(tcpa >= 0.0)
		{
			//Draw CPA line
			bmp_output.draw_line(  cpaRowOwn, cpaColOwn,
								   cpaRowTarget,
								   cpaColTarget,
								   cpaColor   );
		}
	}


	//Text annotations
	for(unsigned int v = 0; v<targetVessels.size(); ++v)
	{
		char txt[1024];

		int vesselRow =  bmp_output.getHeight() / 2 - targetVessels.at(v).getY();
		int vesselCol =  bmp_output.getWidth() / 2  + targetVessels.at(v).getX();


		uint32_t vesselColor = 0x0000FF;

		if((targetVessels.at(v).getTCPA() > 0.0) && (targetVessels.at(v).getTCPA() < 30.0) && (targetVessels.at(v).getCPA() < 50.0))
		{
			vesselColor = 0xFF0000;
		}

		sprintf(txt, "ID: %d", targetVessels.at(v).getID());
		bmp_output.draw_text(txt, vesselRow + 6, vesselCol - 32, vesselColor);
		sprintf(txt, "X: %1.1fm", targetVessels.at(v).getX());
		bmp_output.draw_text(txt, vesselRow + 14, vesselCol - 32, vesselColor);
		sprintf(txt, "Y: %1.1fm", targetVessels.at(v).getY());
		bmp_output.draw_text(txt, vesselRow + 22, vesselCol - 32, vesselColor);
		sprintf(txt, "V: %1.1fm/s", targetVessels.at(v).getSpeed());
		bmp_output.draw_text(txt, vesselRow + 30, vesselCol - 32, vesselColor);
		sprintf(txt, "C: %1.1f", targetVessels.at(v).getCourse());
		bmp_output.draw_text(txt, vesselRow + 38, vesselCol - 32, vesselColor);


		if(targetVessels.at(v).getTCPA() >= 0.0)
		{
			sprintf(txt, "CPA: %1.1fm", targetVessels.at(v).getCPA());
			bmp_output.draw_text(txt, vesselRow + 46, vesselCol - 32, vesselColor);
			sprintf(txt, "TCPA: %1.1fs", targetVessels.at(v).getTCPA());
			bmp_output.draw_text(txt, vesselRow + 54, vesselCol - 32, vesselColor);
		}
		else
		{
			sprintf(txt, "CPA: --");
			bmp_output.draw_text(txt, vesselRow + 46, vesselCol - 32, vesselColor);
			sprintf(txt, "TCPA: --");
			bmp_output.draw_text(txt, vesselRow + 54, vesselCol - 32, vesselColor);
		}
	}


	bmp_output.img_write("radar.bmp");



	char tmpTxt[1024];
	cJSON_AddStringToObject(json_report, "Version", "0.1");
	cJSON_AddStringToObject(json_report, "InputFile", argv[1]);

	sprintf(tmpTxt, "%1.1f", ownVessel.getSpeed());
	cJSON_AddStringToObject(json_report, "OwnVesselSpeed", tmpTxt);


	for(unsigned int v=0; v<targetVessels.size(); ++v)
	{
		sprintf(tmpTxt, "TargetVessel%d",  (v+1));

		cJSON* json_vessel = cJSON_AddObjectToObject(json_report,tmpTxt);

		sprintf(tmpTxt, "%d",  targetVessels.at(v).getID());
		cJSON_AddStringToObject(json_vessel, "ID", tmpTxt);

		sprintf(tmpTxt, "%1.2f",  targetVessels.at(v).getX());
		cJSON_AddStringToObject(json_vessel, "X", tmpTxt);

		sprintf(tmpTxt, "%1.2f",  targetVessels.at(v).getY());
		cJSON_AddStringToObject(json_vessel, "Y", tmpTxt);

		sprintf(tmpTxt, "%1.2f",  targetVessels.at(v).getSpeed());
		cJSON_AddStringToObject(json_vessel, "Speed", tmpTxt);

		sprintf(tmpTxt, "%1.2f",  targetVessels.at(v).getCourse());
		cJSON_AddStringToObject(json_vessel, "Course", tmpTxt);


		if(targetVessels.at(v).getTCPA() >= 0.0)
		{
			sprintf(tmpTxt, "%1.2f",  targetVessels.at(v).getCPA());
			cJSON_AddStringToObject(json_vessel, "CPA", tmpTxt);

			sprintf(tmpTxt, "%1.2f",  targetVessels.at(v).getTCPA());
			cJSON_AddStringToObject(json_vessel, "TCPA", tmpTxt);
		}
		else
		{
			cJSON_AddStringToObject(json_vessel, "CPA", "--");
			cJSON_AddStringToObject(json_vessel, "TCPA", "--");
		}
	}

	FILE* fJSON = fopen("radar.json","w");

	fprintf(fJSON, cJSON_Print(json_report));

	fclose(fJSON);
	cJSON_free(json_report);



	return 0;
}


void inline __attribute__((always_inline)) img_draw_pixel(bmp_img* img, int row, int col, unsigned int rgb)
{
	if(row < 0)
		return;

	if(row >= img->img_header.biHeight)
		return;

	if(col < 0)
		return;

	if(col >= img->img_header.biWidth)
		return;


	bmp_pixel* px = (img->img_pixels[row] + col);

	px->red =   (unsigned char)(rgb >> 16);
	px->green = (unsigned char)(rgb >> 8);
	px->blue =  (unsigned char)(rgb >> 0);
}




void img_draw_digit(bmp_img* img, int digit, int row, int col, unsigned int rgb)
{
	int digitIdx = digit * 8;
	for(int c=0; c<8; ++c)
	{
		for(int r=0; r<8; ++r)
		{
			if( Standard8x8[digitIdx] & (0x01 << r)  )
			{
				img_draw_pixel(img, row + r, col + c, rgb);
			}
		}
		digitIdx++;
	}
}



void BMP_Image::draw_symbol(uint8_t* drawTable, char sym, int row, int col, unsigned int rgb)
{
	for(int r=0; r<8; ++r)
	{
		for(int c=0; c<8; ++c)
		{
			if( *(drawTable + sym*8 + r) & (0x80 >> c)  )
			{
				draw_pixel(row + r, col + c, rgb);
			}
		}
		//digitIdx++;
	}

}



void img_draw_number(bmp_img* img, int number, int row, int col, unsigned int rgb)
{
	char temp[40];

	snprintf(temp, (sizeof(temp)-1), "%d", number);

	for(int i=0; i<(int)strlen(temp); ++i)
	{
		img_draw_digit(img,(temp[i]-48),row, col, rgb);
		col += 8;
	}
}


void BMP_Image::draw_text(const char* txt, int row, int col, unsigned int rgb)
{
	//char temp[40];

	//snprintf(temp, (sizeof(temp)-1), "%d", number);

	int symNum = 0;

	int l = (int)strlen(txt);
	for(int i=0; i<l; ++i)
	{
		draw_symbol((uint8_t*)font8x8_ib8x8u_full, txt[symNum++], row,col,rgb);
		col += 8;
	}
}



void  BMP_Image::draw_line(int row1, int col1, int row2, int col2, unsigned int rgb)
{
	int deltaRow = row2 - row1;
	int deltaCol = col2 - col1;

	int len = (int)(roundf(sqrtf((float)(deltaCol*deltaCol + deltaRow*deltaRow))));

	float stepRow, stepCol;

	float rowFloat, colFloat;
	int rowInt, colInt;

	stepRow = (float)deltaRow / (float)len;
	stepCol = (float)deltaCol / (float)len;

	rowFloat = (float)row1;
	colFloat = (float)col1;

	while(len--)
	{
		rowInt = (int)roundf(rowFloat);
		colInt = (int)roundf(colFloat);

		draw_pixel(rowInt, colInt, rgb);

		rowFloat += stepRow;
		colFloat += stepCol;
	}
}



void  BMP_Image::draw_circle(int row, int col, float radius, unsigned int rgb)
{

	double angle = 0;

	int coordRow = row + (int)(radius * sin(angle));
	int coordCol = col + (int)(radius * cos(angle));;

	while(angle < 360.0)
	{
		angle += 0.5;

		int coordRowNew = row + (int)round(radius * sin( angle * M_PI/180.0 ));
		int coordColNew = col + (int)round(radius * cos( angle * M_PI/180.0 ));;

		draw_line(coordRow, coordCol, coordRowNew, coordColNew, rgb);

		coordRow = coordRowNew;
		coordCol = coordColNew;
	}
}




void inline __attribute__((always_inline)) BMP_Image::draw_pixel(int row, int col, unsigned int rgb)
{
	if(row < 0)
		return;

	if(row >= img.img_header.biHeight)
		return;

	if(col < 0)
		return;

	if(col >= img.img_header.biWidth)
		return;

	bmp_pixel* px = (img.img_pixels[row] + col);

	px->red =   (unsigned char)(rgb >> 16);
	px->green = (unsigned char)(rgb >> 8);
	px->blue =  (unsigned char)(rgb >> 0);
}



void  BMP_Image::img_write(const char* filename)
{
	bmp_img_write(&img, filename);
}
