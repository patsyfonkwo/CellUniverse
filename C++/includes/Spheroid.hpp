#ifndef SPHEROID_HPP
#define SPHEROID_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <array>
#include <algorithm>
#include <string>
#include <unordered_map>
#include "Cell.hpp"
#define M_PI 3.14159265358979323846

class SpheroidParams : public CellParams
{
public:
    float x;
    float y;
    float z;
    std::vector<float> x_vec;
    std::vector<float> y_vec;
    std::vector<float> z_vec;
    float majorRadius;
    float minorRadius;
    float theta_x; // x fixed
    float theta_y; // y fixed
    float theta_z; // z fixed

    SpheroidParams() : CellParams(""), x(0), y(0), z(0), majorRadius(0), minorRadius(0) {}
    SpheroidParams(const std::string &name, float x, float y, float z, float majorRadius, float minorRadius) 
        : CellParams(""), x(x), y(y), z(z), majorRadius(majorRadius), minorRadius(majorRadius) {}
    SpheroidParams(const std::string &name, float x, float y, float z, std::vector<float> _x_vec, std::vector<float> _y_vec, std::vector<float> _z_vec)
        : CellParams(name), x(x), y(y), z(z), x_vec(_x_vec), y_vec(_y_vec), z_vec(_z_vec)
        {
            // assuming this is the x and y vectors are the same and z is the smaller one
            majorRadius = std::sqrt((_x_vec[0]*_x_vec[0])+(_x_vec[1]*_x_vec[1])+(_x_vec[2]*_x_vec[2]));
            minorRadius = std::sqrt((_z_vec[0]*_z_vec[0])+(_z_vec[1]*_z_vec[1])+(_z_vec[2]*_z_vec[2]));
            
            float theta_x = std::atan(z/y); // x fixed
            float theta_y = std::atan(z/x); // y fixed
            float theta_z = std::atan(y/x); // z fixed
        }

    void parseParams(float x_, float y_, float z_, std::vector<float> _x_vec, std::vector<float> _y_vec, std::vector<float> _z_vec)
    {
        x = x_;
        y = y_;
        z = z_;
        x_vec = _x_vec;
        y_vec = _y_vec;
        z_vec = _z_vec;
        majorRadius = std::sqrt((_x_vec[0]*_x_vec[0])+(_x_vec[1]*_x_vec[1])+(_x_vec[2]*_x_vec[2]));
        minorRadius = std::sqrt((_z_vec[0]*_z_vec[0])+(_z_vec[1]*_z_vec[1])+(_z_vec[2]*_z_vec[2]));
    }

    

};

class Spheroid 
{
    private:
        std::vector<std::vector<std::vector<int>>> matrix;
        std::vector<std::vector<std::vector<int>>> rotated_matrix;

        std::string _name;
        std::vector<double> _x_vec;
        std::vector<double> _y_vec;
        std::vector<double> _z_vec;
        cv::Point3f _position;
        double _major_radius;
        double _minor_radius;
        double _rotation;
        double a, b, c;
        bool dormant;


        // double _get_magnitude(std::vector<double> vec);

        bool _is_spheroid (int x, int y, int z, float a, float b, float c);

        std::vector<int> _rotate_point(int x, int y, int z, double theta_x, double theta_y, double theta_z, size_t n);

        std::vector<std::vector<std::vector<int>>> _rotate_matrix(std::vector<std::vector<std::vector<int>>> matrix);
    
    public:
        static SpheroidParams paramClass;
        static SpheroidConfig cellConfig;
        
        Spheroid(const SpheroidParams &init_props);
        
        Spheroid() : _major_radius(0), _minor_radius(0), _rotation(0), dormant(false) {}

        void printCellInfo() {
            std::cout << "Sphere name: " << _name << " x: " << _position.x << " y: " << _position.y << " z: " << _position.z << " majorRadius: " << _major_radius << " minorRadius: " << _minor_radius << " isDormant: " << dormant << std::endl;
        }

        std::vector<double> getShapeAt(double z) const;

        void draw(cv::Mat &image, SimulationConfig simulationConfig, cv::Mat *cellMap = nullptr, float z = 0) const;

        void drawOutline(cv::Mat &image, float color, float z = 0) const;

        [[nodiscard]] Spheroid getPerturbedCell() const;

        Spheroid getParameterizedCell(std::unordered_map<std::string, float> params) const;

        std::tuple<Spheroid, Spheroid, bool> getSplitCells(const std::vector<cv::Mat> &image) const;
        // std::tuple<Sphere, Sphere, bool> Sphere::getSplitCells() const;

        std::vector<std::pair<float, cv::Vec3f>> performPCA(const std::vector<cv::Point3f> &points) const;
        //std::vector<float> performPCA(const std::vector<cv::Point3f> &points) const;

        bool checkConstraints() const;

        CellParams getCellParams() const;

        [[nodiscard]] std::pair<std::vector<float>, std::vector<float>> calculateCorners() const;

        bool checkIfCellsValid(const std::vector<Spheroid> &spheroids)
        {
            return !checkIfCellsOverlap(spheroids);
        }

        std::pair<std::vector<float>, std::vector<float>> calculateMinimumBox(Spheroid &perturbed_cell) const;

        static bool checkIfCellsOverlap(const std::vector<Spheroid> &spheroids);

        float major_magnitude();
        
        float minor_magnitude();

        cv::Point3f get_center() const;

        void print();

        bool paintPixelUpright(int i, int j, int k) const ;

        bool paintPixelRotated(int i, int j, int k) const;

        int get_matrix_size();

       
            

            
        
};

#endif


