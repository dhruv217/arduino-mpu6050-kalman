#include "Wire.h"
#include "I2Cdev.h"
#include "MPU6050.h"

MPU6050 accelgyro;

unsigned long now, lastTime = 0;
float dt;                                   //Differential time

int16_t ax, ay, az, gx, gy, gz;             //Accelerometer gyroscope raw data
float aax=0, aay=0,aaz=0, agx=0, agy=0, agz=0;    //Angle variable
long axo = 0, ayo = 0, azo = 0;             //Accelerometer offset
long gxo = 0, gyo = 0, gzo = 0;             //Gyro offset


float pi = 3.1415926;
float AcceRatio = 16384.0;                  //Accelerometer scale factor
float GyroRatio = 131.0;                    //Gyroscope scale factor

uint8_t n_sample = 8;                       //Accelerometer filter algorithm sampling number
float aaxs[8] = {0}, aays[8] = {0}, aazs[8] = {0};         //x,y-axis sampling queue
long aax_sum, aay_sum,aaz_sum;                      //x,y-axis sampling add

float a_x[10]={0}, a_y[10]={0},a_z[10]={0} ,g_x[10]={0} ,g_y[10]={0},g_z[10]={0}; //Accelerometer covariance calculation queue
float Px=1, Rx, Kx, Sx, Vx, Qx;             //x-axis Calman variables
float Py=1, Ry, Ky, Sy, Vy, Qy;             //y-axis Calman variables
float Pz=1, Rz, Kz, Sz, Vz, Qz;             //z-axis Calman variables

void setup()
{
    Wire.begin();
    Serial.begin(115200);

    accelgyro.initialize();                 //initialization

    unsigned short times = 200;             //The number of samples
    for(int i=0;i<times;i++)
    {
        accelgyro.getMotion6(&ax, &ay, &az, &gx, &gy, &gz); //Read the six-axis original value
        axo += ax; ayo += ay; azo += az;      //Sampling
        gxo += gx; gyo += gy; gzo += gz;
    
    }
    
    axo /= times; ayo /= times; azo /= times; //Calculate accelerometer offset
    gxo /= times; gyo /= times; gzo /= times; //Calculate the gyro offset
}

void loop()
{
    unsigned long now = millis();             //current time(ms)
    dt = (now - lastTime) / 1000.0;           //Differential time(s)
    lastTime = now;                           //Last sampling time(ms)

    accelgyro.getMotion6(&ax, &ay, &az, &gx, &gy, &gz); //Read the six-axis original value

    float accx = ax / AcceRatio;              //x-axis acceleration
    float accy = ay / AcceRatio;              //y-axis acceleration
    float accz = az / AcceRatio;              //z-axis acceleration

    aax = atan(accy / accz) * (-180) / pi;    //The x-axis angle to the z-axis
    aay = atan(accx / accz) * 180 / pi;       //The y-axis angle to the z-axis
    aaz = atan(accz / accy) * 180 / pi;       //The z-axis angle to the y-axis

    aax_sum = 0;                              // Sliding weight filtering algorithm for accelerometer raw data
    aay_sum = 0;
    aaz_sum = 0;
  
    for(int i=1;i<n_sample;i++)
    {
        aaxs[i-1] = aaxs[i];
        aax_sum += aaxs[i] * i;
        aays[i-1] = aays[i];
        aay_sum += aays[i] * i;
        aazs[i-1] = aazs[i];
        aaz_sum += aazs[i] * i;
    
    }
    
    aaxs[n_sample-1] = aax;
    aax_sum += aax * n_sample;
    aax = (aax_sum / (11*n_sample/2.0)) * 9 / 7.0; //Angle AM ​​to 0-90 °
    aays[n_sample-1] = aay;                        //Here we use the experimental method to obtain the appropriate coefficient
    aay_sum += aay * n_sample;                     //This example factor is 9/7
    aay = (aay_sum / (11*n_sample/2.0)) * 9 / 7.0;
    aazs[n_sample-1] = aaz; 
    aaz_sum += aaz * n_sample;
    aaz = (aaz_sum / (11*n_sample/2.0)) * 9 / 7.0;

    float gyrox = - (gx-gxo) / GyroRatio * dt; //x-axis angular velocity
    float gyroy = - (gy-gyo) / GyroRatio * dt; //x-axis angular velocity
    float gyroz = - (gz-gzo) / GyroRatio * dt; //x-axis angular velocity
    agx += gyrox;                             //x-axis angular velocity integral
    agy += gyroy;                             //y-axis angular velocity integral
    agz += gyroz;                             //z-axis angular velocity integral
    
    /* kalman start */
    Sx = 0; Rx = 0;
    Sy = 0; Ry = 0;
    Sz = 0; Rz = 0;
    
    for(int i=1;i<10;i++)
    {                 //The average value of the calculation
        a_x[i-1] = a_x[i];                      //The acceleration average
        Sx += a_x[i];
        a_y[i-1] = a_y[i];
        Sy += a_y[i];
        a_z[i-1] = a_z[i];
        Sz += a_z[i];
    
    }
    
    a_x[9] = aax;
    Sx += aax;
    Sx /= 10;                                 //x-axis acceleration average
    a_y[9] = aay;
    Sy += aay;
    Sy /= 10;                                 //y-axis acceleration average
    a_z[9] = aaz;
    Sz += aaz;
    Sz /= 10;                                 //z-axis acceleration average

    for(int i=0;i<10;i++)
    {
        Rx += sq(a_x[i] - Sx);
        Ry += sq(a_y[i] - Sy);
        Rz += sq(a_z[i] - Sz);
    
    }
    
    Rx = Rx / 9;                              //Get the variance
    Ry = Ry / 9;                        
    Rz = Rz / 9;
  
    Px = Px + 0.0025;                         //0.0025 in the following instructions ...
    Kx = Px / (Px + Rx);                      //Calculate the Kalman gain
    agx = agx + Kx * (aax - agx);             //Gyro angle and accelerometer speed superimposed
    Px = (1 - Kx) * Px;                       //Update p value

    Py = Py + 0.0025;
    Ky = Py / (Py + Ry);
    agy = agy + Ky * (aay - agy); 
    Py = (1 - Ky) * Py;
  
    Pz = Pz + 0.0025;
    Kz = Pz / (Pz + Rz);
    agz = agz + Kz * (aaz - agz); 
    Pz = (1 - Kz) * Pz;

    /* kalman end */

    Serial.print(agx);Serial.print(",");
    Serial.print(agy);Serial.print(",");
    Serial.print(agz);Serial.println();
    
}
