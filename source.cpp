
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <fstream>
#include <string>
#include <string.h>
#include <ctime>
#include <direct.h>
#include <random>

//TO DO: redo all the commentaries in english


using namespace std;

// параметры симул€ции
double desiredTemperature;
int boxSize, amountOfParticles;
string filename;

// random
random_device rd;
mt19937 gen(rd());

// физические параметры
const double sigma = 1;
const double initialMass = 1;
const double epsilon = 1;
// sigma = 3.4 * 10^-10 м; initialMass (дл€ аргона) = 6.63 * 10^-32 кг; epsilon (дл€ аргона) = 1.65 * 10^-21 ƒж
// 1 sqrt(initialMass * epsilon / sigma / sigma) = 2.15 пс (единица времени в программе)
const double sigmaReal = 3.4 * pow(10, -10);
const double initialMassReal = 6.63 * pow(10, -26);
const double epsilonReal = 1.65 * pow(10, -21);
const double k_b = 1.38 * pow(10, -23);


// параметры программы
const double timeStep = 0.001; // 2.15 фс
const double cutOffDistance = sigma * 2.5;
const double zeroPotentialEnergy = 4 * epsilon * (pow(sigma / cutOffDistance, 12) - pow(sigma / cutOffDistance, 6));
const double collisionFrequency = 1.5;
const int endTime = 100;
const double relaxationTimeShare = 0.4;
const int pressureArraySize = (int)floor(endTime / timeStep * (1 - relaxationTimeShare));
double pressure[pressureArraySize];

// класс частиц
class Particle {
public:
    double Fx, Fy, Fz;
    double vx, vy, vz;
    double x, y, z;
    double xp, yp, zp;
    double xt, yt, zt;
    double mass;
    double difx, dify, difz;

    Particle() {
        mass = 0;
        Fx = Fy = Fz = vx = vy = vz = x = y = z = 0;
    }

    Particle(double coordinates[3], double speed[3], double aMass, double time) {
        x = coordinates[0];
        y = coordinates[1];
        z = coordinates[2];
        vx = speed[0];
        vy = speed[1];
        vz = speed[2];
        xp = coordinates[0] - vx * time;
        yp = coordinates[1] - vy * time;
        zp = coordinates[2] - vz * time;
        mass = aMass;
        difx = dify = difz = 0;
    }

    void move(double time) {
        xt = x;
        yt = y;
        zt = z;
        x = 2 * x - xp + Fx * pow(time, 2);
        y = 2 * y - yp + Fy * pow(time, 2);
        z = 2 * z - zp + Fz * pow(time, 2);
        vx = (x - xp) / 2 / time;
        vy = (y - yp) / 2 / time;
        vz = (z - zp) / 2 / time;
        xp = xt;
        yp = yt;
        zp = zt;
        Fx = Fy = Fz = 0;
        difx += vx * time;
        dify += vy * time;
        difz += vz * time;
    }



};

// удержание координат частицы в пределах коробки
Particle normalizeCoordinates(Particle p, double time) {
    double arr[3] = {p.x, p.y, p.z};
    double arrp[3] = {p.xp, p.yp, p.zp};
    for (int i = 0; i < 3; i++) {
        while (arr[i] >= boxSize) {
            arr[i] -= boxSize;
            arrp[i] -= boxSize;
        }
        while (arr[i] < 0) {
            arr[i] += boxSize;
            arrp[i] += boxSize;
        }
    }
    p.x = arr[0];
    p.y = arr[1];
    p.z = arr[2];
    p.xp = arrp[0];
    p.yp = arrp[1];
    p.zp = arrp[2];
    return p;
}

// скал€рное произведение
double scalar(double x[3], double y[3]) {
    double result = 0;
    for (int i = 0; i < 3; i++) {
        result += x[i] * y[i];
    }
    return result;
}


// измеренные величины
double potentialEnergy = 0;

void getFilename() {
    cout << "Enter name of the save:" << endl;
    cin >> filename;
    const char* filepath = ("saves/" + filename).c_str();
    const char* savespath = "saves";
    mkdir(savespath);
    mkdir(filepath);
}

void getParameters() {
    getFilename();
    cout << "Enter desired temperature:" << endl;
    cin >> desiredTemperature;
    cout << "Enter size of a box in sigma units: " << endl;
    cin >> boxSize;
    cout << "Enter amount of particles: " << endl;
    cin >> amountOfParticles;
}




double* calculateMinimalDistance(double relativeCoordinates[3]) {
    for (int i = 0; i < 3; i++) {
        if (abs(relativeCoordinates[i]) > boxSize/2) {
            relativeCoordinates[i] = (boxSize-abs(relativeCoordinates[i])) * (-relativeCoordinates[i]) / abs(relativeCoordinates[i]);
        }
    }
    return relativeCoordinates;

}

Particle* createParticles(const int amountOfParticles, Particle* particles) {
    bool hasEnoughDistance = true;
    double minimalDistance = 10 * sigma;
    for (int i = 0; i < amountOfParticles; i++) {
        double aCoordinates[3];
        aCoordinates[0] = ((double)rand() / RAND_MAX * boxSize);
        aCoordinates[1] = (double)rand() / RAND_MAX * boxSize;
        aCoordinates[2] = (double)rand() / RAND_MAX * boxSize;
        double aSpeed[3];
        aSpeed[0] = (-1 + (double)rand() / RAND_MAX * 2);
        aSpeed[1] = (-1 + (double)rand() / RAND_MAX * 2);
        aSpeed[2] = (-1 + (double)rand() / RAND_MAX * 2);
        particles[i] = Particle(aCoordinates, aSpeed, initialMass, timeStep);
        minimalDistance = 10 * sigma;
        for (int j = 0; j < i; j++) {
            double relativeCoordinates[3] = {particles[i].x-particles[j].x, particles[i].y-particles[j].y, particles[i].z-particles[j].z};
            calculateMinimalDistance(relativeCoordinates);
            double relativeDistance = scalar(relativeCoordinates, relativeCoordinates);
            minimalDistance = min(minimalDistance, relativeDistance);
            };
        while (minimalDistance < 1.3 * sigma) {
            minimalDistance = 10 * sigma;
            aCoordinates[0] = ((double)rand() / RAND_MAX * boxSize);
            aCoordinates[1] = (double)rand() / RAND_MAX * boxSize;
            aCoordinates[2] = (double)rand() / RAND_MAX * boxSize;
            particles[i] = Particle(aCoordinates, aSpeed, initialMass, timeStep);
            for (int j = 0; j < i; j++) {
                double relativeCoordinates[3] = {particles[i].x-particles[j].x, particles[i].y-particles[j].y, particles[i].z-particles[j].z};
                calculateMinimalDistance(relativeCoordinates);
                double relativeDistance = scalar(relativeCoordinates, relativeCoordinates);
                minimalDistance = min(minimalDistance, relativeDistance);
            };

        };

    };
    return particles;
}



Particle* createParticles(Particle* particles) {
    int t = 3;
    assert(amountOfParticles == t*t*t);
    double x, y, z;
    double aCoordinates[3];
    for (int i = 0; i < t*t*t; i++) {
        x = (double)boxSize / t * ((i / (t* t)) % t);
        y = (double)boxSize / t * ((i / t) % t);
        z = (double)boxSize / t * (i % t);
        double aSpeed[3];
        aSpeed[0] = (-5 + (double)rand() / RAND_MAX * 10);
        aSpeed[1] = (-5 + (double)rand() / RAND_MAX * 10);
        aSpeed[2] = (-5 + (double)rand() / RAND_MAX * 10);
        aCoordinates[0] = x;
        aCoordinates[1] = y;
        aCoordinates[2] = z;
        particles[i] = Particle(aCoordinates, aSpeed, 1.0, timeStep);
    };
    return particles;
}

Particle* create2Particles(Particle* particles) {
    assert(amountOfParticles == 2);
    double x, y, z;
    double aCoordinates[3];
    for (int i = 0; i < 2; i++) {
        x = 4 + i * 5;
        y = 1;
        z = 1;
        double aSpeed[3];
        aSpeed[0] = 2;
        aSpeed[1] = 0;
        aSpeed[2] = 0;
        aCoordinates[0] = x;
        aCoordinates[1] = y;
        aCoordinates[2] = z;
        particles[i] = Particle(aCoordinates, aSpeed, 1.0, timeStep);
    };
    return particles;
}


double* calculateForce(Particle P1, Particle P2, double force[3]) {
    double Fx, Fy, Fz;
    Fx = Fy = Fz = 0.0;
    force[0] = force[1] = force[2] = 0.0;
    double dx = P1.x - P2.x;
    double dy = P1.y - P2.y;
    double dz = P1.z - P2.z;
    double relativeCoordinates[3] = { dx, dy, dz };
    calculateMinimalDistance(relativeCoordinates);
    double relativeDistance = scalar(relativeCoordinates, relativeCoordinates);
    if (relativeDistance > cutOffDistance * cutOffDistance) { return force; }
    double forceAbsolute = epsilon * (48 * pow(sigma, 12) / pow(relativeDistance, 7) - 24 * pow(sigma, 6) / pow(relativeDistance, 4));
    potentialEnergy += 4 * epsilon * (pow(sigma, 12) / pow(relativeDistance, 6) - pow(sigma, 6) / pow(relativeDistance, 3)) - zeroPotentialEnergy;
    force[0] = forceAbsolute * relativeCoordinates[0];
    force[1] = forceAbsolute * relativeCoordinates[1];
    force[2] = forceAbsolute * relativeCoordinates[2];
    assert(sqrt(pow(force[0], 2) + pow(force[1], 2) + pow(force[2], 2)) < 1000000);
    return force;
}

Particle thermostat(Particle p, double time) {
    //Implementation of Andersen thermostat
    double probability = (double) rand() / RAND_MAX;
    if (probability >= time * collisionFrequency) { return p; };
    std::normal_distribution<double> v_dist(0, sqrt(k_b * desiredTemperature / initialMassReal));
    double vx = v_dist(gen) / sqrt(epsilonReal / initialMassReal);
    double vy = v_dist(gen) / sqrt(epsilonReal / initialMassReal);
    double vz = v_dist(gen) / sqrt(epsilonReal / initialMassReal);
    p.vx = vx;
    p.vy = vy;
    p.vz = vz;
    p.xp = p.x - vx * time;
    p.yp = p.y - vy * time;
    p.zp = p.z - vz * time;
    return p;
}

double* impulse(Particle* p, double mv[3]) {
    mv[0] = 0;
    mv[1] = 0;
    mv[2] = 0;
    for (int i = 0; i < amountOfParticles; i++) {
        mv[0] += p[i].vx;
        mv[1] += p[i].vy;
        mv[2] += p[i].vz;
    }
    return mv;
}

double getKineticEnergy(Particle* particles) {
    double K = 0;
    for (int i = 0; i < amountOfParticles; i++) {
        K += pow(particles[i].vx, 2) / 2;
        K += pow(particles[i].vy, 2) / 2;
        K += pow(particles[i].vz, 2) / 2;
    }
    return K * epsilonReal;
}

void calculateImmediatePressure(Particle* p, int time) {
    double virial = 0;
    double V = pow(boxSize * sigmaReal, 3);
    //double pi = acos(-1);
    //double p_tail = 16 / 3 * pi * pow((initialMassReal * amountOfParticles / V), 2) * (2 / 3 * pow((1 / cutOffDistance * sigmaReal),9) - pow((1 / cutOffDistance * sigmaReal),3));
    double K = getKineticEnergy(p);
    for (int i = 0; i < amountOfParticles; i++) {
        double force[3] = {p[i].Fx, p[i].Fy, p[i].Fz};
        double coords[3] = {p[i].x, p[i].y, p[i].z};
        virial += scalar(force, coords);
    }
    virial *= epsilonReal;
    /*cout << virial << endl;
    cout << K * 2 << endl;
    cout << (2 * K + virial) / 3 / V << endl;
    cout << (2 * K) / 3 / V << endl;*/
    //cout << p_tail << endl;
    pressure[time] = (2 * K + virial) / 3 / V;
}

void giveSpeed(Particle* p, double mvInit[3], double timeStep) {
    std::normal_distribution<double> v_dist(0, sqrt(k_b * desiredTemperature / initialMassReal));
    double vx, vy, vz;
    vx = vy = vz = 0;
    for (int i = 0; i < amountOfParticles; i++) {
        p[i].vx = v_dist(gen) / sqrt(epsilonReal / initialMassReal);;
        p[i].vy = v_dist(gen) / sqrt(epsilonReal / initialMassReal);;
        p[i].vz = v_dist(gen) / sqrt(epsilonReal / initialMassReal);;
    }
    impulse(p, mvInit);
    double temperature = getKineticEnergy(p) * 2 / 3 / amountOfParticles / k_b;
    for (int i = 0; i < amountOfParticles; i++) {
        p[i].vx -= mvInit[0] / amountOfParticles;
        p[i].vy -= mvInit[1] / amountOfParticles;
        p[i].vz -= mvInit[2] / amountOfParticles;

        /*p[i].vx *= sqrt(desiredTemperature / temperature);
        p[i].vx *= sqrt(desiredTemperature / temperature);
        p[i].vx *= sqrt(desiredTemperature / temperature);*/

        p[i].xp = p[i].x - p[i].vx * timeStep;
        p[i].yp = p[i].y - p[i].vy * timeStep;
        p[i].zp = p[i].z - p[i].vz * timeStep;

        p[i].difx = p[i].dify = p[i].difz = 0;
    }
    impulse(p, mvInit);
    cout << mvInit[0] << " " << mvInit[1] << " " << mvInit[2] << endl;
}

void checkExistence(string filepath) {
    ifstream fin;
    fin.open(filepath);
    if (fin.is_open()) {
        return;
    }
    else {
        fin.close();
        ofstream fout;
        fout.open(filepath);
        fout << amountOfParticles << endl;
        fout << endl;
        fout.close();
    };
}

void writeEnergy(int timer, Particle* particles) {
    if (timer % 10 != 0) { return; };
    ofstream fout;
    string filepath = "./saves/" + filename + "/_energy.txt";
    fout.open(filepath, ios_base::app);
    assert(fout.is_open());
    double K = getKineticEnergy(particles);
    double U = potentialEnergy * epsilonReal;
    fout << timer << " " << K << " " << U << endl;
    fout.close();
}

void writeCoordinates(int timer, Particle* particles) {
    if (timer % 50 != 0) { return; };
    ofstream fout;
    string filepath = "./saves/" + filename + "/" + to_string(timer) + ".txt";
    checkExistence(filepath);
    fout.open(filepath);
    assert(fout.is_open());
    fout << amountOfParticles << endl;
    fout << endl;
    for (int i = 0; i < amountOfParticles; i++) {
        fout << "Ar " << particles[i].x << " " << particles[i].y << " " << particles[i].z << endl;
    }
    fout.close();
}

void checkDiffusionFileExistence(string filepath) {
    ifstream fin;
    fin.open(filepath);
    if (fin.is_open()) {
        return;
    }
    else {
        fin.close();
        ofstream fout;
        fout.open(filepath);
        fout << sigmaReal << " " << sigmaReal * sqrt(initialMassReal / epsilonReal) * timeStep << endl;
        fout.close();
    };
}

void writeDiffusion(int timer, Particle* p) {
    if (timer % 100 != 0) { return; };
    ofstream fout;
    string filepath = "./saves/" + filename + "/_diffusion.txt";
    checkDiffusionFileExistence(filepath);
    fout.open(filepath, ios_base::app);
    assert(fout.is_open());
    fout << timer << " ";
    double diffusion = 0;
    for (int i = 0; i < amountOfParticles; i ++) {
        diffusion += (pow(p[i].difx, 2) + pow(p[i].dify, 2) + pow(p[i].difz, 2));
    }
    diffusion /= amountOfParticles;
    fout << diffusion << endl;
    fout.close();
}

void writeSpeed(int timer, Particle* particles, double mvInit[3]) {
    if (timer % 50 != 0) { return; };
    ofstream fout;
    string filepath = "./saves/" + filename + "/_speed.txt";
    fout.open(filepath, ios_base::trunc);
    assert(fout.is_open());
    double k = sqrt(epsilonReal / initialMassReal);
    for (int i = 0; i < amountOfParticles; i++) {
        fout << to_string(particles[i].vx * k) + " " + to_string(particles[i].vy * k) + " " + to_string(particles[i].vz * k) << endl;
    }
    fout.close();

    string absFilepath = "./saves/" + filename + "/_speedAbs.txt";
    fout.open(absFilepath, ios_base::app);
    assert(fout.is_open());
    double absSpeed;
    for (int i = 0; i < amountOfParticles; i++) {
        absSpeed = sqrt(pow(particles[i].vx*k, 2) + pow(particles[i].vy*k, 2) + pow(particles[i].vz*k, 2));
        fout << absSpeed;
        if (i < amountOfParticles-1) {
            fout << " ";
        }
    }
    fout << endl;
    fout.close();
}

void writeFinalInformation(Particle* p, int timer, double mvInit[3]) {
    double finalPressure = 0;
    double mvFin[3];
    impulse(p, mvFin);
    for (int i = 0; i < pressureArraySize; i++) {
        finalPressure += pressure[i] / pressureArraySize;
    }
    ofstream fout;
    string filepath = "./saves/" + filename + "/_info.txt";
    fout.open(filepath);
    assert(fout.is_open());
    fout << "Save name: " << filename << endl;
    fout << endl;
    fout << "Sigma unit in meters: " << sigmaReal << endl;
    fout << "Particle mass in kilograms: " << initialMassReal << endl;
    fout << "Epsilon unit in Joule: " << epsilonReal << endl;
    fout << "Time unit in seconds: " << sigmaReal * sqrt(initialMassReal / epsilonReal) << endl;
    fout << endl;
    fout << "Amount of particles: " << amountOfParticles << endl;
    fout << "Length of a cube box in sigma units: " << boxSize << endl;
    fout << "Desired temperature in K: " << desiredTemperature << endl;
    fout << endl;
    fout << "Pressure in Pa: " << finalPressure << endl;
    fout << "Temperature as mean kinetic energy (in K): " << getKineticEnergy(p) * 2 / 3 / amountOfParticles / k_b << endl;
    fout << "Initial system impulse in reduced units: " << mvInit[0] << " " << mvInit[1] << " " << mvInit[2] << endl;
    fout << "Final system impulse in reduced units: " << mvFin[0] << " " << mvFin[1] << " " << mvFin[2] << endl;
    fout << endl;
    fout << "Total time steps in simulation:" << endTime / timeStep << endl;
    fout << "Total time in simulation time units:" << endTime << endl;
    fout << "Length of one timestep in seconds:" << sigmaReal * sqrt(initialMassReal / epsilonReal) * timeStep << endl;
    fout << "Total time passed:" << sigmaReal * sqrt(initialMassReal / epsilonReal)* endTime << endl;
    fout << "Timesteps spent on relaxation: " << endTime * relaxationTimeShare / timeStep << endl;
    fout << "Time spent on relaxation: " << sigmaReal * sqrt(initialMassReal / epsilonReal) * endTime * relaxationTimeShare << endl;
    fout << "Average amount of collisons per time unit with heat bath in relaxation period:" << collisionFrequency << endl;
    fout.close();
}

int main()
{
    srand(time(0));
    setlocale(LC_ALL, "ru");
    int timer = 0;

    double temperature = 0;
    double mvInit[3];
    double mvFin[3];

    getParameters();

    Particle* particles = new Particle[amountOfParticles];
    createParticles(amountOfParticles, particles);

    bool hasNulledImpulse = false;

    double vx, vy, vz;
    vx = vy = vz = 0.0;
    for (int i = 0; i < amountOfParticles; i++) {
        vx += particles[i].vx;
        vy += particles[i].vy;
        vz += particles[i].vz;
    }
    for (int i = 0; i < amountOfParticles; i++) {
        particles[i].vx -= vx / amountOfParticles;
        particles[i].vy -= vy / amountOfParticles;
        particles[i].vz -= vz / amountOfParticles;
    }

    while (true) {

        potentialEnergy = 0;

        for (int i = 0; i < amountOfParticles; i++) {
            for (int j = i + 1; j < amountOfParticles; j++) {
                double dForce[3];
                calculateForce(particles[i], particles[j], dForce);
                particles[i].Fx += dForce[0];
                particles[i].Fy += dForce[1];
                particles[i].Fz += dForce[2];
                particles[j].Fx -= dForce[0];
                particles[j].Fy -= dForce[1];
                particles[j].Fz -= dForce[2];
            }
        }

        if (!(timer * timeStep < endTime * relaxationTimeShare)) {
            calculateImmediatePressure(particles, (timer - int(endTime * relaxationTimeShare / timeStep)));
        }

        for (int i = 0; i < amountOfParticles; i++) {
            particles[i].move(timeStep);
            particles[i] = normalizeCoordinates(particles[i], timeStep);
            if (timer * timeStep < endTime * relaxationTimeShare) {
                    particles[i] = thermostat(particles[i], timeStep);
            }
        }

        timer++;
        if (timer * timeStep > endTime) {
            writeFinalInformation(particles, timer, mvInit);
            return 0;
        };

        if (!(timer * timeStep < endTime * relaxationTimeShare)) {
            if (!hasNulledImpulse) {
                cout << "Nullifying impulse" << endl;
                giveSpeed(particles, mvInit, timeStep);
                hasNulledImpulse = true;
                impulse(particles, mvFin);

            }
            writeEnergy(timer, particles);
            writeSpeed(timer, particles, mvInit);
            writeDiffusion(timer, particles);
        }

        writeCoordinates(timer, particles);

        cout << timer << endl;
    }

    return 0;
}
