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
const double collisionFrequency = 0.5;

// класс частиц
class Particle {
public:
    double Fx, Fy, Fz;
    double vx, vy, vz;
    double x, y, z;
    double xp, yp, zp;
    double xt, yt, zt;
    double mass;

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
        Fx = 0;
        Fy = 0;
        Fz = 0;
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

void thermostat(Particle p, double time) {
    //Implementation of Andersen thermostat
    double probability = (double) rand() / RAND_MAX;
    if (probability >= time * collisionFrequency) { return; };
    std::normal_distribution<double> v_dist(0, sqrt(k_b * desiredTemperature / initialMassReal));
    double vx = v_dist(gen) / sqrt(epsilonReal / initialMassReal);
    double vy = v_dist(gen) / sqrt(epsilonReal / initialMassReal);
    double vz = v_dist(gen) / sqrt(epsilonReal / initialMassReal);
    p.vx = vx;
    p.vy = vy;
    p.vz = vz;
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

void writeSpeed(int timer, Particle* particles) {
    if (timer % 100 != 0) { return; };
    ofstream fout;
    string filepath = "./saves/" + filename + "/_speed.txt";
    fout.open(filepath, ios_base::trunc);
    assert(fout.is_open());
    for (int i = 0; i < amountOfParticles; i++) {
        fout << to_string(particles[i].vx) + " " + to_string(particles[i].vy) + " " + to_string(particles[i].vz) << endl;
    }
    fout.close();
}

int main()
{

    srand(time(0));
    setlocale(LC_ALL, "ru");
    int timer = 0;

    double temperature = 0;

    getParameters();

    Particle* particles = new Particle[amountOfParticles];
    createParticles(amountOfParticles, particles);

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

        for (int i = 0; i < amountOfParticles; i++) {
            particles[i].move(timeStep);
            particles[i] = normalizeCoordinates(particles[i], timeStep);
            thermostat(particles[i], timeStep);
        }

        writeEnergy(timer, particles);
        writeSpeed(timer, particles);
        writeCoordinates(timer, particles);

        timer++;
        if (timer > 50000) { return 0; };
        cout << timer << endl;
    }

    return 0;
}
