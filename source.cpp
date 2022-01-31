#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

//TO DO: redo all the commentaries in english


using namespace std;

// ��������� ���������
double initialTemperature;
int boxSize, amountOfParticles;
string filename;

// ���������� ���������
const double sigma = 1;
const double initialMass = 1;
const double epsilon = 1;
// sigma = 3.4 * 10^-10 �; initialMass (��� ������) = 6.63 * 10^-32 ��; epsilon (��� ������) = 1.65 * 10^-21 ��
// 1 sqrt(initialMass * epsilon / sigma / sigma) = 2.15 �� (������� ������� � ���������)


// ��������� ���������
const double timeStep = 0.001; // 2.15 ��
const double cutOffDistance = sigma * 4;
const double zeroPotentialEnergy = 4 * epsilon * (pow(sigma / cutOffDistance, 12) - pow(sigma / cutOffDistance, 6));

// ��������� ��������� ������� � �������� �������
double* normalizeCoordinates(double coordinates[3]) {
    for (int i = 0; i < sizeof(coordinates) / sizeof(coordinates[0]); i++) {
        while (coordinates[i] >= boxSize) {
            coordinates[i] -= boxSize;
        }
        while (coordinates[i] < 0) {
            coordinates[i] += boxSize;
        }
    }
    return coordinates;
}

// ��������� ������������
double scalar(double x[3], double y[3]) {
    double result = 0;
    for (int i = 0; i < 3; i++) {
        result += x[i] * y[i];
    }
    return result;
}


// ���������� ��������
double potentialEnergy = 0;

void getFilename() {}

void getParameters() {
    getFilename();
    cout << "������� ��������� �����������:" << endl;
    cin >> initialTemperature;
    cout << "������� ������ �������: " << endl;
    cin >> boxSize;
    cout << "������� ���������� ������: " << endl;
    cin >> amountOfParticles;
}


// ����� ������
class Particle {
public:
    double Fx, Fy, Fz;
    double vx, vy, vz;
    double x, y, z;
    double mass;

    Particle() {
        mass = 0;
        Fx = Fy = Fz = vx = vy = vz = x = y = z = 0;
    }

    Particle(double coordinates[3], double speed[3], double aMass) {
        x = coordinates[0];
        y = coordinates[1];
        z = coordinates[2];
        vx = speed[0];
        vy = speed[1];
        vz = speed[2];
        mass = aMass;
    }

    void move(double time) {
        vx += Fx * time;
        vy += Fy * time;
        vz += Fz * time;
        x += vx * time;
        y += vy * time;
        z += vz * time;
    }



};

Particle* createParticles(const int amountOfParticles) {
    Particle* particles = new Particle[amountOfParticles];
    for (int i = 0; i < amountOfParticles; i++) {
        double aCoordinates[3];
        aCoordinates[0] = ((double)rand() / RAND_MAX * boxSize);
        aCoordinates[1] = (double)rand() / RAND_MAX * boxSize;
        aCoordinates[2] = (double)rand() / RAND_MAX * boxSize;
        double aSpeed[3];
        aSpeed[0] = (-1 + (double)rand() / RAND_MAX * 2);
        aSpeed[1] = (-1 + (double)rand() / RAND_MAX * 2);
        aSpeed[2] = (-1 + (double)rand() / RAND_MAX * 2);
        particles[i] = Particle(aCoordinates, aSpeed, initialMass);
    }
    return particles;
}

double* calculateMinimalDistance(double relativeCoordinates[3]) {
    for (int i = 0; i < 3; i++) {
        if (abs(relativeCoordinates[i]) > abs(relativeCoordinates[i] - boxSize)) {
            relativeCoordinates[i] += (-boxSize) * relativeCoordinates[i] / abs(relativeCoordinates[i]);
        }
    }
    return relativeCoordinates;

}

double* calculateForce(Particle P1, Particle P2, double force[3]) {
    double Fx, Fy, Fz;
    Fx = Fy = Fz = 0;
    double dx = P1.x - P2.x;
    double dy = P1.y - P2.y;
    double dz = P1.z - P2.z;
    double relativeCoordinates[3] = { dx, dy, dz };
    calculateMinimalDistance(relativeCoordinates);
    double relativeDistance = scalar(relativeCoordinates, relativeCoordinates);
    if (relativeDistance > cutOffDistance * cutOffDistance) { return force; }
    double forceAbsolute = epsilon * (48 * pow(sigma, 12) / pow(relativeDistance, 7) - 24 * pow(sigma, 6) / pow(relativeDistance, 4));
    potentialEnergy += 4 * epsilon * (pow(sigma, 12) / pow(relativeDistance, 6) - pow(sigma, 6) / pow(relativeDistance, 3));
    force[0] = forceAbsolute * dx;
    force[1] = forceAbsolute * dy;
    force[2] = forceAbsolute * dz;
    return force;
}

void writeEnergy(double timer, Particle* particles) {

}

void writeSpeed(double timer, Particle* particles) {

}

int main()
{
    setlocale(LC_ALL, "ru");
    getParameters();
    Particle* particles = createParticles(amountOfParticles);
    int timer = 0;
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
        }
        writeEnergy(timer, particles);
        writeSpeed(timer, particles);
        timer++;
        cout << timer << endl;
    }

    return 0;
}
