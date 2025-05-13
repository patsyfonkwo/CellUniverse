#include "../../includes/Spheroid.hpp"
#include <random>
// #include <iostream>

SpheroidParams Spheroid::paramClass = SpheroidParams();
SpheroidConfig Spheroid::cellConfig = SpheroidConfig();



static double _get_magnitude(std::vector<double> vec){
    return std::sqrt((vec[0]*vec[0])+(vec[1]*vec[1])+(vec[2]*vec[2]));
}

bool Spheroid::_is_spheroid (int x, int y, int z, float a, float b, float c){
    // builds matrix; returns true if matrix coordinates should be considered in the spheroid
    // if (((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c) < 1.1) && ((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c) > 0.9))
    if (((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c) < 1.1))
        return true;
    return false;
}

std::vector<int> Spheroid::_rotate_point(int x, int y, int z, double theta_x, double theta_y, double theta_z, size_t n){

    // will return std::vector<int> of rotated points
    std::vector<std::vector<double>> Rx = {{1, 0, 0},
                                            {0, std::cos(theta_x), -std::sin(theta_x)},
                                            {0, std::sin(theta_x), std::cos(theta_x)}};

    std::vector<std::vector<double>> Ry = {{std::cos(theta_y), 0, std::sin(theta_y)}, 
                                            {0, 1, 0},
                                            {-std::sin(theta_y), 0, std::cos(theta_y)}};

    std::vector<std::vector<double>> Rz = {{std::cos(theta_z), -std::sin(theta_z), 0},
                                            {std::sin(theta_z), std::cos(theta_z), 0},
                                            {0, 0, 1}};

    
    // std::cout << "OG x : " << x << std::endl;
    // std::cout << "OG y : " << y << std::endl;
    // std::cout << "OG z : " << z << std::endl;

    // std::cout << "theta_x : " << theta_x << std::endl;
    // std::cout << "cos(theta_x) : " << std::cos(theta_x) << std::endl;

    // Rotate around x
    double new_x = (double)x;
    double new_y = (double)y * std::cos(theta_x) + (double)z * std::sin(theta_x);
    double new_z = (double)z * std::cos(theta_x) - (double)y * std::sin(theta_x);

    // std::cout << "x after x : " << new_x << std::endl;
    // std::cout << "y after x : " << new_y << std::endl;
    // std::cout << "z after x : " << new_z << std::endl;

    // Rotate around y
    new_x = new_x * std::cos(theta_y) - new_z * std::sin(theta_y);
    new_y = new_y;
    new_z = new_x * std::sin(theta_y) + new_z * std::cos(theta_y);

    // std::cout << "x after y : " << new_x << std::endl;
    // std::cout << "y after y : " << new_y << std::endl;
    // std::cout << "z after y : " << new_z << std::endl;

    // Rotate around z
    new_x = new_x * std::cos(theta_z) + new_y * std::sin(theta_z);
    new_y = new_y * std::cos(theta_z) - new_x * std::sin(theta_z);;
    new_z = new_z;

    // std::cout << "x after z : " << new_x << std::endl;
    // std::cout << "y after z : " << new_y << std::endl;
    // std::cout << "z after z : " << new_z << std::endl;

    std::vector<int> new_point;
    if ((new_x < 0) || (new_x >= n) || (new_y < 0) || (new_y >= n) || (new_z < 0) || (new_z >= n)){
        new_point = {-1, -1, -1}; // point rotated off matrix
        // std::cout << "Off the matrix!" << std::endl;
        // std::cout << "N : " << n << std::endl;
    }
    else {
        new_point = {(int)std::floor(new_x), (int)std::floor(new_y), (int)std::floor(new_z)};
        // std::cout << "On the matrix!" << std::endl;
        // std::cout << "x after z : " << new_x << std::endl;
        // std::cout << "y after z : " << new_y << std::endl;
        // std::cout << "z after z : " << new_z << std::endl;
    }    
    
    
    return new_point;

}

std::vector<std::vector<std::vector<int>>> Spheroid::_rotate_matrix(std::vector<std::vector<std::vector<int>>> matrix){
    // new dimensions derived from angles of x, y, and z
    // Angles about axes: theta_x, theta_y, theta_z
    float x = this->_major_radius; // assuming x and y
    float y = this->_major_radius; // vectors are majors (oblate)
    float z = this->_minor_radius;

    float theta_x = std::atan(z/y) ; // x fixed
    float theta_y = std::atan(z/x); // y fixed
    float theta_z = std::atan(y/x); // z fixed


    // int matrix_size = std::ceil(_major_radius * 2);
    std::vector<std::vector<std::vector<int>>> rotated_matrix((int)matrix.size()*2, std::vector<std::vector<int>>((int)matrix.size()*2, std::vector<int>((int)matrix.size()*2, 0)));
    
    // offset values to make up for (0, 0, 0) indexing
    

    int offset = rotated_matrix.size()/2;
    // rotate each point on matrix onto rotated_matrix by all three angles
    for (size_t i = 0; i < matrix.size(); ++i){
        for (size_t j = 0; j < matrix.size(); ++j){
            for (size_t k = 0; k < matrix.size(); ++k){
                // treat point (i, j, k) like a vector being rotated about center
                std::vector<int> new_point = _rotate_point(i+offset, j+offset, k+offset, theta_x, theta_y, theta_z, matrix.size());
                // if ((new_point[0]+offset>=0) && (new_point[0]+offset<rotated_matrix.size()) && (new_point[1]+offset>=0) && (new_point[1]+offset<rotated_matrix.size()) && (new_point[2]+offset>=0) && (new_point[2]+offset<rotated_matrix.size()))
                //     rotated_matrix[new_point[0]+offset][new_point[1]+offset][new_point[2]+offset] = matrix[i][j][k];
                if ((new_point[0]-offset>=0) && (new_point[0]-offset<rotated_matrix.size()) && (new_point[1]-offset>=0) && (new_point[1]-offset<rotated_matrix.size()) && (new_point[2]-offset>=0) && (new_point[2]-offset<rotated_matrix.size()))
                    rotated_matrix[new_point[0]-offset][new_point[1]-offset][new_point[2]-offset] = matrix[i][j][k];
                // else 
                    // std::cout << "not here" << std::endl;
            }
        }
    }
    
    return rotated_matrix;
}

Spheroid::Spheroid(const SpheroidParams &init_props) 
: _name(init_props.name), _position{init_props.x, init_props.y, init_props.z},
          _major_radius(init_props.majorRadius), _minor_radius(init_props.majorRadius), _rotation(0), dormant(false)
{
    // std::cout << "NAME : " << _name << std::endl;
    // std::cout << "XPOS : " << _position.x << std::endl;
    // std::cout << "YPOS : " << _position.y << std::endl;
    // std::cout << "ZPOS : " << _position.z << std::endl;
    // std::cout << "Major : " << _major_radius << std::endl;
    // std::cout << "Minor : " << _minor_radius << std::endl;
    // std::cout << "Init Props Major : " << init_props.majorRadius << std::endl;
    // std::cout << "Init Props Minor : " << init_props.minorRadius << std::endl;
    this->a = this->_major_radius; // major-radius could also be average of both vector magnitudess
    this->b = this->a; // _get_magnitude(_y_vec); // process directly
    this->c = this->_minor_radius;

    // use radii to build formula for 
    
    // assuming oblate

    // initialize matrix 
    std::vector<float> xyz = {std::ceil(2 * a), std::ceil(2 * b), std::ceil(2 * c)};

    int matrix_size = *std::max_element(xyz.begin(), xyz.end()) + 10 ; // 10 pixel padding
    // std::cout << "matrix_size : " << matrix_size << std::endl; 
    // std::cout << "_major_radius : " << _major_radius << std::endl; 
    // std::cout << "_minor_radius : " << _minor_radius << std::endl; 

    std::vector<std::vector<std::vector<int>>> matrix(matrix_size, std::vector<std::vector<int>>(matrix_size, std::vector<int>(matrix_size, 0)));
    for (int i=0; i < matrix_size; ++i){
        for (int j=0; j < matrix_size; ++j){
            for (int k=0; k < matrix_size; ++k){
                if (_is_spheroid(i-a, j-b, k-c, a, b, c)){
                    matrix[i][j][k] = 1;
                    // std::cout << "(" << i << ", " << j << ", " << j << ")" << std::endl;
                }
                // std::cout << matrix[i][j][k];
            }
            // std::cout << std::endl;
        }
    }
    this->matrix = matrix;

// print rotated matrix
    // this->rotated_matrix = _rotate_matrix(matrix);
    // for (int i=0; i < rotated_matrix.size(); ++i){
    //     for (int j=0; j < rotated_matrix.size(); ++j){
    //         for (int k=0; k < rotated_matrix.size(); ++k){
    //             std::cout << rotated_matrix[i][j][k];
    //         }
    //         std::cout << std::endl;
    //     }
    // }
    // int xSize = std::ceil(2 * a);
    // int ySize = std::ceil(2 * b);
    // int zSize = std::ceil(2 * c);
    // std::vector<std::vector<std::vector<int>>> matrix(xSize, std::vector<std::vector<int>>(ySize, std::vector<int>(zSize, 0)));
    // for (int i=0; i < xSize; ++i){
    //     for (int j=0; j < ySize; ++j){
    //         for (int k=0; k < zSize; ++k){
    //             if (_is_spheroid(i-a, j-b, k-c, a, b, c)){
    //                 matrix[i][j][k] = 1;
    //             }
    //             std::cout << matrix[i][j][k];
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    // std::cout << std::endl;
}

float Spheroid::major_magnitude(){
    return _get_magnitude(_x_vec); // semi-major radius
}

float Spheroid::minor_magnitude(){
    return _get_magnitude(_z_vec); // semi-minor radius
}

cv::Point3f Spheroid::get_center() const {
    return _position; // x, y, z position
}

void Spheroid::print(){
    for (int i=0; i < this->matrix.size(); ++i){
        for (int j=0; j < this->matrix.size(); ++j){
            for (int k=0; k < this->matrix.size(); ++k){
                std::cout << matrix[i][j][k];
            }
            std::cout << std::endl;
        }
    }
}

bool Spheroid::paintPixelUpright(int i, int j, int k) const {
    if ((k >= matrix[0][0].size()) || (j >= matrix[0].size()) || (i >= matrix.size()))
        return false;
    return (matrix[i][j][k] == 1);
}

bool Spheroid::paintPixelRotated(int i, int j, int k) const {
    if ((k >= rotated_matrix[0][0].size()) || (j >= rotated_matrix[0].size()) || (i >= rotated_matrix.size()))
        return false;
    return (rotated_matrix[i][j][k] == 1);
}

int Spheroid::get_matrix_size(){
    return matrix.size();
}

// bool is_ellipse(double x, double y, double z, std::vector<std::vector<double>> R){
//     if (((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c) < 1.1) && ((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c) > 0.9))
//         return true;
//     return false;
// }

std::vector<double> Spheroid::getShapeAt(double z) const
{
    // need to find axes and radii using z value
    std::vector<double> vec = {_major_radius, _minor_radius, paramClass.x, paramClass.y};
    return vec;
}

void Spheroid::draw(cv::Mat &image, SimulationConfig simulationConfig, cv::Mat *cellMap, float z) const{
    if (dormant)
    {
        return;
    }

    std::vector<double> currentShape = getShapeAt(z); // major, minor, x, and y

    if ((currentShape[0] <= 0) || (currentShape[1] <= 0))
    {
        std::cout << "cell not found" << std::endl; // this keeps on executing (draws when commented out)
        return; // cell doesn't exist or unseen by algorithm
    }

    if ((currentShape[2] < 0) || (currentShape[3] < 0) || (currentShape[2] >= 517) || (currentShape[3] >= 422))
    {
        // std::cout << "cell not in frame" << std::endl; // this doesn't execute
        return; // x and y values not on the piture
    }

    // paint rotated pixels onto the z slice 1 by 1
    int xpos = get_center().x;
    int ypos = get_center().y;
    int zpos = get_center().z;

    // std::cout << "NAME : " << _name << std::endl;
    // std::cout << "XPOS : " << xpos << std::endl;
    // std::cout << "YPOS : " << ypos << std::endl;
    // std::cout << "ZPOS : " << zpos << std::endl;

    // int matrix_offset = get_matrix_size() / 2;
    
    for (int i=0; (i+xpos < image.cols); ++i){
        for (int j=0; (j+ypos < image.rows); ++j){
            if (paintPixelUpright(i-xpos, j-ypos, z - zpos)){
                std::cout << "should paint : " << "(" << j+ypos << ", " << i+xpos << ", " << z << ")" << std::endl;
                // image.at<cv::Vec3b>(i+x_offset, j+y_offset) = cv::Vec3b(255, 255, 255);
                image.at<float>(j+ypos, i+xpos) = simulationConfig.cell_color; 
            }
            // image.at<float>(j+y_offset, i+x_offset) = 0.0f; // simulationConfig.cell_color;
        }
    }
    // if (image.empty())
    //     std::cout << "empty" << std::endl;
    // std::cout << "Image Columns: " << image.cols << std::endl;
    // std::cout << "Image Rows: " << image.rows << std::endl;
    // std::cout << "Image Size: " << image.size() << std::endl;

    // for (int i=0; (i < image.rows); ++i){
    //     for (int j=0; (j < image.cols); ++j){
    //         image.at<float>(i, j) = 0.0f;
    //     }
    // }

    // // float background_color = simulationConfig.background_color;
    // // float cell_color = simulationConfig.cell_color;

    // // cv::Point center(_position.x, _position.y);
    // // cv::circle(image, center, static_cast<int>(currentRadius), cv::Scalar(cell_color), -1);
    // // Define ellipse parameters

    // cv::Point center(currentShape[2], currentShape[2]);
    // cv::Size axes(currentRadii.first, currentRadii.second); // major and minor axes
    // double angle = 45.0;

    // double startAngle = 0.0;
    // double endAngle = 360.0;

    // cv::Scalar color(255, 255, 255); // White color
    // int thickness = 3;
    // int lineType = cv::LINE_AA;

    // cv::ellipse(image, center, axes, angle, startAngle, endAngle, color, thickness, lineType);
}

void Spheroid::drawOutline(cv::Mat &image, float color, float z) const {
    int x_offset = get_center().x;
    int y_offset = get_center().y;
    int z_offset = get_center().z;

    cv::Vec3b color_vec(255, 255, 255); // White color
    // int matrix_offset = get_matrix_size() / 2;
    
    // for (int i=0; (i+x_offset < image.size().width); ++i){
    //     for (int j=0; (j+y_offset < image.size().height); ++j){
    //         if (paintPixelUpright(i-x_offset, j-y_offset, z - z_offset)){
    //             // image.at<cv::Vec3b>(i + x_offset, j + y_offset) = cv::Vec3b(0, 0, 0);
    //             // image.at<float>(i + x_offset, j + y_offset) = 0;
    //             // cv::Vec3b* pixel = this->tiff_slices[k+z_offset].ptr<cv::Vec3b>(j+y_offset) + i+x_offset;
    //             // *pixel = cv::Vec3b(color);
    //         }
    //         image.at<cv::Vec3b>(i + x_offset, j + y_offset) = color_vec;
    //     }
    // }
    // for (int i=0; (i < 517); ++i){
    //     for (int j=0; (j < 422); ++j){
    //         image.at<cv::Vec3b>(i, j) = color_vec; // simulationConfig.cell_color;
    //     }
    // }
}

[[nodiscard]] Spheroid Spheroid::getPerturbedCell() const {
    SpheroidParams spheroidParams(
        _name,
        // FIXME: we should choose only ONE of these, uniformly at random, to perturb in each iteration.
        _position.x + cellConfig.x.getPerturbOffset(),
        _position.y + cellConfig.y.getPerturbOffset(),
        _position.z + cellConfig.z.getPerturbOffset(),
        _major_radius + cellConfig.majorRadius.getPerturbOffset(),
        _minor_radius + cellConfig.minorRadius.getPerturbOffset());
    return Spheroid(spheroidParams);
}

Spheroid Spheroid::getParameterizedCell(std::unordered_map<std::string, float> params) const {
    float xOffset = params["x"];
    float yOffset = params["y"];
    float zOffset = params["z"];
    float majorRadiusOffset = params["majorRadius"];
    float minorRadiusOffset = params["minorRadius"];

    if (params.empty())
    {
        xOffset = Spheroid::cellConfig.x.getPerturbOffset();
        yOffset = Spheroid::cellConfig.y.getPerturbOffset();
        zOffset = Spheroid::cellConfig.z.getPerturbOffset();
        majorRadiusOffset = Spheroid::cellConfig.majorRadius.getPerturbOffset();
        minorRadiusOffset = Spheroid::cellConfig.minorRadius.getPerturbOffset();
    }

    float newMajorRadius = fmin(fmax(Spheroid::cellConfig.minMajorRadius, _major_radius + majorRadiusOffset), Spheroid::cellConfig.maxMajorRadius);
    float newMinorRadius = fmin(fmax(Spheroid::cellConfig.minMinorRadius, _minor_radius + minorRadiusOffset), Spheroid::cellConfig.maxMinorRadius);
    SpheroidParams spheroidParams(
        _name,
        _position.x + xOffset,
        _position.y + yOffset,
        _position.z + zOffset,
        newMajorRadius,
        newMinorRadius);
    return Spheroid(spheroidParams);
}

std::tuple<Spheroid, Spheroid, bool> Spheroid::getSplitCells(const std::vector<cv::Mat> &image) const {
    // remove calback function
    // make z scale same as x and y scale
    // 16 z slices per x-y
    // interpolate (make up) z slices convert 400x400x30 to 400x400x400
    // using brightness values interpolate
    // Step 1: Get the bounding box using calculateCorners
    // Step 1: Get the bounding box using calculateCorners
    // print out eigenvectors for just x,y,z directions
    // for the z vector, we can scale it instead of relying on interpolation
    auto [min_corner, max_corner] = calculateCorners();

    int minX = std::max(0, static_cast<int>(min_corner[0]));
    int maxX = std::min(static_cast<int>(image[0].cols), static_cast<int>(max_corner[0]));
    int minY = std::max(0, static_cast<int>(min_corner[1]));
    int maxY = std::min(static_cast<int>(image[0].rows), static_cast<int>(max_corner[1]));
    int minZ = std::max(0, static_cast<int>(min_corner[2]));
    int maxZ = std::min(static_cast<int>(image.size()), static_cast<int>(max_corner[2]));

    // Step 2: Construct 3D points for PCA
    std::vector<cv::Point3f> points;
    for (int z = minZ; z <= maxZ; ++z) {
        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                points.emplace_back(x, y, z);
                
            }
        }
    }
    if (!points.empty()) {
        // Perform PCA to get eigenvalues
        std::vector<std::pair<float, cv::Vec3f>> eigenvalues = performPCA(points);

        // Print the top 3 eigenvalues
        std::cout << "Eigenvalues with most variance:" << std::endl;
        for (size_t i = 0; i < eigenvalues.size(); ++i) {
            std::cout << "Eigenvalue " << i + 1 << ": " << eigenvalues[i].first << " " << eigenvalues[i].second << std::endl;
        }
    } else {
        std::cout << "No points found for PCA." << std::endl;
    }

    double theta = ((double)rand() / RAND_MAX) * 2 * M_PI;
    double phi = ((double)rand() / RAND_MAX) * M_PI;

    cv::Point3f split_axis(
        sin(phi) * cos(theta),
        sin(phi) * sin(theta),
        cos(phi));

    cv::Point3f offset = split_axis * (_major_radius / 2.0);
    cv::Point3f new_position1 = _position + offset;
    cv::Point3f new_position2 = _position - offset;

    double halfMajorRadius = _major_radius / 2.0;
    double halfMinorRadius = _minor_radius / 2.0;

    Spheroid cell1(SpheroidParams(_name + "0", new_position1.x, new_position1.y, new_position1.z, halfMajorRadius, halfMinorRadius));
    Spheroid cell2(SpheroidParams(_name + "1", new_position2.x, new_position2.y, new_position2.z, halfMajorRadius, halfMinorRadius));

    bool constraints = cell1.checkConstraints() && cell2.checkConstraints();

    // Cell cell1Base = static_cast<Cell>(cell1);
    // Cell cell2Base = static_cast<Cell>(cell2);//convert Sphere to Cell

    return std::make_tuple(Spheroid(cell1), Spheroid(cell2), constraints);
}
// std::tuple<Sphere, Sphere, bool> Sphere::getSplitCells() const;

std::vector<std::pair<float, cv::Vec3f>> Spheroid::performPCA(const std::vector<cv::Point3f> &points) const {
    if (points.empty())
    {
        throw std::invalid_argument("No points provided for PCA.");
    }

    // Create a matrix from the points
    cv::Mat data(points.size(), 3, CV_32F);
    for (size_t i = 0; i < points.size(); ++i)
    {
        data.at<float>(i, 0) = points[i].x;
        data.at<float>(i, 1) = points[i].y;
        data.at<float>(i, 2) = points[i].z;
    }

    // Perform PCA
    cv::PCA pca(data, cv::Mat(), cv::PCA::DATA_AS_ROW);

    // Extract the eigenvalues and eigenvectors
    cv::Mat eigenvalues = pca.eigenvalues;
    cv::Mat eigenvectors = pca.eigenvectors;

    // Prepare the result
    std::vector<std::pair<float, cv::Vec3f>> eigenPairs;
    for (int i = 0; i < std::min(3, eigenvalues.rows); ++i)
    {
        float eigenvalue = eigenvalues.at<float>(i);
        cv::Vec3f eigenvector(
            eigenvectors.at<float>(i, 0),
            eigenvectors.at<float>(i, 1),
            eigenvectors.at<float>(i, 2)
        );
        eigenPairs.emplace_back(eigenvalue, eigenvector);
    }

    return eigenPairs;
}
//std::vector<float> performPCA(const std::vector<cv::Point3f> &points) const;

bool Spheroid::checkConstraints() const {
    SpheroidConfig config;
    return (config.minMajorRadius <= _major_radius) && (_major_radius <= config.maxMajorRadius) && (config.minMinorRadius <= _minor_radius) && (_minor_radius <= config.maxMinorRadius);
}

CellParams Spheroid::getCellParams() const {
    return SpheroidParams(_name, _position.x, _position.y, _position.z, _major_radius, _minor_radius);
}

[[nodiscard]] std::pair<std::vector<float>, std::vector<float>> Spheroid::calculateCorners() const {
    std::vector<float> min_corner = {static_cast<float>(_position.x) - static_cast<float>(_major_radius),
                                     static_cast<float>(_position.y) - static_cast<float>(_major_radius),
                                     static_cast<float>(_position.z) - static_cast<float>(_major_radius)};

    std::vector<float> max_corner = {static_cast<float>(_position.x) + static_cast<float>(_major_radius),
                                     static_cast<float>(_position.y) + static_cast<float>(_major_radius),
                                     static_cast<float>(_position.z) + static_cast<float>(_major_radius)};

    return std::make_pair(min_corner, max_corner);
}

std::pair<std::vector<float>, std::vector<float>> Spheroid::calculateMinimumBox(Spheroid &perturbed_cell) const {
    auto [cell1_min_corner, cell1_max_corner] = calculateCorners();
    auto [cell2_min_corner, cell2_max_corner] = perturbed_cell.calculateCorners();

    std::vector<float> min_corner, max_corner;
    for (int i = 0; i < 3; ++i)
    {
        min_corner.push_back(std::min(cell1_min_corner[i], cell2_min_corner[i]));
        max_corner.push_back(std::max(cell1_max_corner[i], cell2_max_corner[i]));
    }
    return std::make_pair(min_corner, max_corner);
}

bool Spheroid::checkIfCellsOverlap(const std::vector<Spheroid> &spheroids) {
    std::vector<std::vector<float>> positions;
    std::vector<std::pair<float, float>> radii;

    for (const auto &cell : spheroids)
    {
        positions.push_back({cell._position.x, cell._position.y, cell._position.z});
        radii.push_back({cell._major_radius * 0.95, cell._minor_radius * 0.95});
    }

    std::vector<std::vector<float>> distances;
    for (const auto &position1 : positions)
    {
        std::vector<float> distance_row;
        for (const auto &position2 : positions)
        {
            float distance = 0.0f;
            for (int i = 0; i < 3; ++i)
            {
                distance += pow(position1[i] - position2[i], 2);
            }
            distance = sqrt(distance);
            distance_row.push_back(distance);
        }
        distances.push_back(distance_row);
    }

    std::vector<std::vector<std::pair<float, float>>> radii_sums;
    for (const auto &radius1 : radii)
    {
        std::vector<std::pair<float, float>> radii_row;
        for (const auto &radius2 : radii)
        {
            radii_row.push_back({radius1.first + radius2.first, radius1.second + radius2.second});
        }
        radii_sums.push_back(radii_row);
    }

    bool overlap = false;
    for (std::size_t i = 0; i < spheroids.size(); ++i)
    {
        for (std::size_t j = 0; j < spheroids.size(); ++j)
        {
            if (i != j && (distances[i][j] < radii_sums[i][j].first) && (distances[i][j] < radii_sums[i][j].second))
            {
                overlap = true;
                break;
            }
        }
        if (overlap)
        {
            break;
        }
    }

    return overlap;
}