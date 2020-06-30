#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <typeinfo>

const double EPSILON = 0.0001;

enum Sign{
    NEGATIVE,
    POSITIVE
};

enum Turn{
    LEFT,
    RIGHT,
    COLLINEAR
};

bool areEqual(const double& first, const double& second) {
    return abs(first - second) <= EPSILON;
}



/////////////////POINT STRUCTURE/////////////////
struct Point {
public:
    Point() : x(0.0), y(0.0) {};
    Point(double x, double y) : x(x), y(y) {};
    Point(const Point&) = default;

    bool operator == (const Point&) const;
    bool operator != (const Point&) const;

    Point& operator +=(const Point&);
    Point& operator -=(const Point&);

    Point operator+(const Point&) const;
    Point operator-(const Point&) const;

    double x;
    double y;
};

bool Point::operator==(const Point& other) const{
    return areEqual(x, other.x) && areEqual(y, other.y);
}

bool Point::operator!=(const Point& other) const{
    return !(*this == other);
}

Point& operator *=(Point& lhs, const double& rhs) {
    lhs.x *= rhs;
    lhs.y *= rhs;
    return lhs;
}

Point& operator *=(const double& lhs, Point& rhs) {
    rhs.x *= lhs;
    rhs.y *= lhs;
    return rhs;
}

Point& operator /=(Point& lhs, const double& rhs) {
    lhs.x /= rhs;
    lhs.y /= rhs;
    return lhs;
}

Point operator *(const Point& lhs, const double& rhs) {
    Point result = lhs;
    return result *= rhs;
}

Point operator *(const double& lhs, const Point& rhs) {
    Point result = rhs;
    return result *= lhs;
}

Point operator /(const Point& lhs, const double& rhs) {
    Point result = lhs;
    return result /= rhs;
}


Point& Point::operator+=(const Point& other) {
    x += other.x;
    y += other.y;
    return *this;
}

Point& Point::operator-=(const Point& other) {
    x -= other.x;
    y -= other.y;
    return *this;
}

Point Point::operator+(const Point& other) const {
    Point result = *this;
    return result += other;
}

Point Point::operator-(const Point& other) const {
    Point result = *this;
    return result -= other;
}

double getSqrtDistance(const Point& first, const Point& second) {
    return sqrt((first.x - second.x) * (first.x - second.x) + (first.y - second.y) * (first.y - second.y));
}

double getDistance(const Point& first, const Point& second) {
    return (first.x - second.x) * (first.x - second.x) + (first.y - second.y) * (first.y - second.y);
}

double crossProduct(const Point& first, const Point& second) {
    return first.x * second.y - second.x * first.y;
}

Turn orientation(const Point& first, const Point& second, const Point& third) {
    double turn = crossProduct(second - first, third - first);
    if (turn > EPSILON) {
        return LEFT;
    }
    if (turn < -EPSILON) {
        return RIGHT;
    }
    return COLLINEAR;
}



/////////////////LINE CLASS/////////////////
class Line {
public:
    Line() : a(0), b(0), c(0) {};
    Line(const Point& first, const Point& second) : a(first.y - second.y), b(second.x - first.x),
                                                    c(first.x * second.y - second.x * first.y) {toNormalState();};
    Line(const double& slope, const double& shift) : a(slope), b (-1.0), c(shift) {};
    Line(const Point& point, const double& slope) : a(slope), b(-1.0), c(point.y - slope*point.x) {toNormalState();};
    Line (const double& a, const double& b, const double& c) : a(a), b(b), c(c) {toNormalState();};
    Line(const Line&) = default;

    bool operator ==(const Line&) const;
    bool operator != (const Line&) const;

    void toNormalState();
    double distanceToPoint(const Point&) const;

    double a, b, c;
};

bool Line::operator ==(const Line& other) const {
    return areEqual(a * other.b, b * other.a) && areEqual(b * other.c, c * other.b);
}

bool Line::operator !=(const Line& other) const {
    return !(*this == other);
}

double crossProduct(const Line& first, const Line& second) {
    return -first.b * second.a + first.a * second.b;
}

void Line::toNormalState() {
    if (!areEqual(b, 0.)) {
        a /= -b;
        b /= -b;
        c /= -b;
    }
}

double Line::distanceToPoint(const Point& point) const{
    return std::abs(point.x * a + point.y * b + c) / sqrt(a * a + b * b);
}




/////////////////SHAPE CLASS/////////////////
class Shape{
public:
    virtual ~Shape() = default;

    virtual double perimeter() const = 0;
    virtual double area() const = 0;

    virtual bool operator ==(const Shape&) const = 0;
    virtual bool operator !=(const Shape&) const = 0;

    virtual bool isCongruentTo(const Shape&) const = 0;
    virtual bool isSimilarTo(const Shape&) const = 0;
    virtual bool containsPoint(const Point&) const = 0;

    virtual void rotate(const Point&, const double&) = 0;
    virtual void reflex(const Point&) = 0;
    virtual void reflex(const Line&) = 0;
    virtual void scale(const Point&, const double&) = 0;

};




/////////////////ELLIPSE CLASS/////////////////
class Ellipse : public Shape {
public:
    Ellipse() : firstFocus(Point(0., 0.)), secondFocus(Point(0., 0.)), distanceToFocuses(0.),
                ellipseCenter(){};
    Ellipse(const Point&, const Point&, const double&);
    Ellipse(const Ellipse&) = default;

    double perimeter() const override;
    double area() const override;

    bool operator==(const Shape&) const override;
    bool operator!=(const Shape&) const override;

    bool isCongruentTo(const Shape&) const override;
    bool isSimilarTo(const Shape&) const override;
    bool containsPoint(const Point&) const override;

    void rotate(const Point&, const double&) override;
    void reflex(const Point&) override;
    void reflex(const Line&) override;
    void scale(const Point&, const double&) override;

    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    Point center() const;

protected:
    Point firstFocus;
    Point secondFocus;
    double distanceToFocuses;
    Point ellipseCenter;

    double majorAxis() const;
    double minorAxis() const;
};

Ellipse::Ellipse(const Point& a, const Point& b, const double& d) {
    firstFocus = a;
    secondFocus = b;
    distanceToFocuses = d;
    ellipseCenter = (a + b) * 0.5;
    std::cerr << firstFocus.x << " " << firstFocus.y << " "<<secondFocus.x << " " << secondFocus.y <<"\n";
    std::cerr << ellipseCenter.x <<" "<<ellipseCenter.y <<"\n";
    std::cerr << distanceToFocuses << "\n";
}

double Ellipse::perimeter() const {
    double a = majorAxis();
    double b = minorAxis();
    return M_PI * (3. * (a + b) - sqrt(3. * a * a + 10. * a * b + 3. * b * b));
}

double Ellipse::area() const {
    double a = majorAxis();
    double b = minorAxis();
    return std::fabs(M_PI * a * b);
}

bool Ellipse::operator==(const Shape& other) const {
    try {
        const Ellipse& another = dynamic_cast<const Ellipse&>(other);
        return (areEqual(distanceToFocuses, another.distanceToFocuses) &&
               ((firstFocus == another.firstFocus && secondFocus == another.secondFocus) ||
               (firstFocus == another.secondFocus && secondFocus == another.firstFocus)));
    } catch (std::bad_cast& error) {
        return false;
    }
}

bool Ellipse::operator!=(const Shape &other) const {
    return !(*this == other);
}

bool Ellipse::isCongruentTo(const Shape& other) const {
    try {

        const Ellipse& another = dynamic_cast<const Ellipse&>(other);

        return areEqual(distanceToFocuses, another.distanceToFocuses) &&
               areEqual(getDistance(firstFocus, secondFocus), getDistance(another.firstFocus, another.secondFocus));
    } catch (std::bad_cast& error) {
        return false;
    }
}

bool Ellipse::isSimilarTo(const Shape& other) const {
    try {
        const Ellipse& another = dynamic_cast<const Ellipse&>(other);
        return areEqual(getSqrtDistance(firstFocus, secondFocus) / getSqrtDistance(another.firstFocus, another.secondFocus),
                        distanceToFocuses / another.distanceToFocuses);
    } catch (std::bad_cast& error) {
        return false;
    }
}

bool Ellipse::containsPoint(const Point& point) const {
    return getSqrtDistance(point, firstFocus) + getSqrtDistance(point, secondFocus) <= distanceToFocuses;
}

void Ellipse::rotate(const Point& center, const double& angle) {
    double newAngle = (angle * M_PI) / 180.0;
    double sinA = sin(newAngle);
    double cosA = cos(newAngle);
    double newfX = (firstFocus.x - center.x) * cosA - (firstFocus.y - center.y) * sinA + center.x;
    double newsY = (secondFocus.x - center.x) * sinA + (secondFocus.y - center.y) * cosA + center.y;
    firstFocus.x = newfX;
    secondFocus.y = newsY;
    ellipseCenter = (firstFocus + secondFocus) * 0.5;
}

void Ellipse::reflex(const Point& center) {
    firstFocus = center + center - firstFocus;
    secondFocus = center + center - secondFocus;
    ellipseCenter = (firstFocus + secondFocus) * 0.5;
}

void Ellipse::reflex(const Line& axis) {
    Point firstPoint(axis.a, axis.b);
    double distToFirstFocus = axis.distanceToPoint(firstFocus);
    double distToSecondFocus = axis.distanceToPoint(secondFocus);
    Point firstShift = firstPoint * (distToFirstFocus / (sqrt(firstPoint.x * firstPoint.x + firstPoint.y * firstPoint.y)));
    Point secondShift = firstPoint * (distToSecondFocus / (sqrt(firstPoint.x * firstPoint.x + firstPoint.y * firstPoint.y)));

    firstFocus += areEqual(axis.distanceToPoint(firstFocus + firstShift * 2), distToFirstFocus) ? firstShift * 2. :
                                                                                                             firstShift * -2.;

    secondFocus += areEqual(axis.distanceToPoint(secondFocus + secondShift * 2), distToSecondFocus) ? secondShift * 2. :
                                                                                                                 secondShift * -2.;
    ellipseCenter = (firstFocus + secondFocus) * 0.5;
}

void Ellipse::scale(const Point& center, const double& coefficient) {
    firstFocus = (firstFocus - center) * coefficient + center;
    secondFocus = (secondFocus - center) * coefficient + center;;
    distanceToFocuses = distanceToFocuses * std::fabs(coefficient);
    ellipseCenter = (firstFocus + secondFocus) * 0.5;
}

std::pair<Point, Point> Ellipse::focuses() const {
    return std::make_pair(firstFocus, secondFocus);
}

std::pair<Line, Line> Ellipse::directrices() const {
    return std::make_pair(Line(eccentricity(), 0.0, -majorAxis()), Line(eccentricity(), 0.0, majorAxis()));
}

double Ellipse::eccentricity() const {
    return getSqrtDistance(firstFocus, secondFocus) / (2. * majorAxis());
}

Point Ellipse::center() const {
    return ellipseCenter;
}

double Ellipse::majorAxis() const {
    return distanceToFocuses * 0.5;
}

double Ellipse::minorAxis() const {
    return sqrt(majorAxis() * majorAxis() - getDistance(firstFocus, secondFocus) / 4.);
}



/////////////////CIRCLE CLASS/////////////////
class Circle : public Ellipse {
public:
    Circle(const Point& circleCenter, const double& circleRadius);

    double perimeter() const override;
    double radius() const;

    void scale(const Point&, const double& coefficient) override;

protected:
    double circleRadius;
};

Circle::Circle(const Point& center, const double& radius) {
    firstFocus = center;
    secondFocus = center;
    ellipseCenter = center;
    circleRadius = radius;
    distanceToFocuses = 2. * radius;
    std::cerr << firstFocus.x << " " << firstFocus.y << " "<<secondFocus.x << " " << secondFocus.y <<"\n";
    std::cerr << distanceToFocuses << "\n";
}

double Circle::perimeter() const {
    return M_PI * 2. * circleRadius;
}

double Circle::radius() const {
    return circleRadius;
}

void Circle::scale(const Point& center, const double& coefficient) {
    firstFocus = (firstFocus - center) * coefficient + center;
    secondFocus = (secondFocus - center) * coefficient + center;;
    distanceToFocuses = distanceToFocuses * std::fabs(coefficient);
    circleRadius = circleRadius * std::fabs(coefficient);
    ellipseCenter = (firstFocus + secondFocus) * 0.5;
}




/////////////////POLYGON CLASS/////////////////
class Polygon : public Shape {
public:
    Polygon() : vertices(std::vector<Point>(0)), verticesNum(0) {};
    Polygon(const std::vector<Point>& vertices) : vertices(vertices), verticesNum(vertices.size()) {};
    template<typename... Args>
    Polygon(Args...);
    Polygon(const Polygon& other) = default;

    int verticesCount() const;
    std::vector<Point> getVertices() const;
    bool isConvex() const;

    double perimeter() const override;
    double area() const  override;

    bool operator ==(const Shape&) const override;
    bool operator !=(const Shape&) const override;

    bool isCongruentTo(const Shape&) const override;
    bool isSimilarTo(const Shape&) const override;
    bool containsPoint(const Point&) const override;

    void rotate(const Point&, const double&) override;
    void reflex(const Point&) override;
    void reflex(const Line&) override;
    void scale(const Point&, const double&) override;

protected:
    std::vector<Point> vertices;
    int verticesNum;

};

template <typename... Args>
Polygon::Polygon(Args... args) {
    for (const Point& point : std::initializer_list<Point>{args...}) {
        vertices.push_back(point);
    }
    verticesNum = vertices.size();
}

int Polygon::verticesCount() const {
    return verticesNum;
}

std::vector<Point> Polygon::getVertices() const {
    return vertices;
}

bool Polygon::isConvex() const {
    Sign orderSign = crossProduct(Line(vertices[0], vertices[1]),
                                  Line(vertices[0], vertices[2])) >= 0. ? POSITIVE : NEGATIVE;
    int size = vertices.size();
    for (int i = 1; i < size - 2; ++i) {
        double crossP = crossProduct(Line(vertices[i], vertices[i + 1]),
                                           Line(vertices[i], vertices[i + 2]));
        Sign crossProductSign = crossP >= 0. ? POSITIVE : NEGATIVE;
        if (crossProductSign != orderSign) {
            return false;
        }
    }
    return true;
}

double Polygon::perimeter() const {
    double perimeter = 0.;
    int size = vertices.size();
    for (int i = 0; i < size; ++i) {
        perimeter += getSqrtDistance(vertices[i], vertices[(i + 1) % size]);
    }
    return perimeter;
}

double Polygon::area() const {
    double area = 0.;
    int size = vertices.size();
    for (int i = 0; i < size; ++i) {
        area += vertices[i].x * vertices[(i + 1) % size].y;
        area -= vertices[i].y * vertices[(i + 1) % size].x;
    }
    return fabs(area) / 2.;
}

bool Polygon::operator==(const Shape &other) const {
    try {
        const Polygon& another = dynamic_cast<const Polygon&>(other);
        if (vertices.size() != another.vertices.size()) {
            return false;
        }
        int size = vertices.size();
        int firstEqual = 0;
        for (int i = 0; i < size; ++i) {
            if (another.vertices[i] == vertices[0]) {
                firstEqual = i;
                break;
            }
        }
        bool rightOrder = true;
        bool leftOrder = true;
        for (int i = 0; i < size; ++i) {
            if (vertices[i] != another.vertices[(firstEqual + i) % size]) {
                rightOrder = false;
                break;
            }
        }
        for (int i = 0; i < size; ++i) {
            if (vertices[i] != another.vertices[(size + firstEqual - i) % size]) {
                leftOrder = false;
                break;
            }
        }
        return leftOrder || rightOrder;
    } catch (std::bad_cast& error) {
        return false;
    }
}

bool Polygon::operator!=(const Shape &other) const {
    return !(*this == other);
}

bool Polygon::isCongruentTo(const Shape &other) const {
    return isSimilarTo(other) && areEqual(area(), other.area());
}

bool Polygon::isSimilarTo(const Shape &other) const {
    try {

        const Polygon& another = dynamic_cast<const Polygon&>(other);

        if (vertices.size() != another.vertices.size()) {
            return false;
        }
        double k = perimeter() / another.perimeter();
        int size = vertices.size();
        int firstEqual = 0;
        for (int i = 0; i < size; ++i) {
            double anotherSide = getSqrtDistance(another.vertices[i], another.vertices[(i + 1) % size]);
            double side = getSqrtDistance(vertices[0], vertices[1]);
            if (areEqual(side / anotherSide, k)) {
                firstEqual = i;
                break;
            }
        }

        bool leftOrder = true;
        bool rightOrder = true;
        for (int i = 0; i < size; ++i) {
            double anotherSide = getSqrtDistance(another.vertices[firstEqual + i],
                                             another.vertices[(firstEqual + i + 1) % size]);
            double side = getSqrtDistance(vertices[i], vertices[(i + 1) % size]);
            if (!areEqual(side / anotherSide, k)) {
                rightOrder = false;
                break;
            }
        }
        for (int i = 0; i < size; ++i) {
            double anotherSide = getSqrtDistance(another.vertices[(size + firstEqual - i + 1)% size],
                                             another.vertices[(size + firstEqual - i) % size]);
            double side = getSqrtDistance(vertices[i], vertices[(i + 1) % size]);
            if (!areEqual(side / anotherSide, k)) {
                leftOrder = false;
                break;
            }
        }
        return leftOrder || rightOrder;
    } catch (std::bad_cast& error) {
        return false;
    }

}

bool Polygon::containsPoint(const Point& point) const {
    int intersectonsNum = 0;
    for (int i = 0; i < verticesNum; ++i) {
        Point start = vertices[i];
        Point end = vertices[(i + 1) % verticesNum];

        if (start.y == end.y) {
            continue;
        }
        if (start.y > end.y) {
            std::swap(start, end);
        }
        if (point.y == std::max(start.y, end.y) && point.x < std::min(start.x, end.x)) {
            ++intersectonsNum;
            continue;
        }
        if (point.y == std::min(start.y, end.y)) {
            continue;
        }

        Turn orient = orientation(start, end, point);
        if (orient == COLLINEAR) {
            return true;
        }
        if (point.y >= std::min(start.y, end.y) && point.y < std::max(start.y, end.y) && orient == LEFT) {
            ++intersectonsNum;
            continue;
        }

    }
    return intersectonsNum % 2;
}

void Polygon::rotate(const Point& center, const double& angle) {
    double newAngle = (angle * M_PI) / 180.0;
    double cosA = cos(newAngle);
    double sinA = sin(newAngle);
    for (Point& vertex: vertices) {
        double newX = (vertex.x - center.x) * cosA - (vertex.y - center.y) * sinA + center.x;
        double newY = (vertex.x - center.x) * sinA + (vertex.y - center.y) * cosA + center.y;
        vertex = Point(newX, newY);
    }

}

void Polygon::reflex(const Point& center) {
    for (Point& vertex : vertices) {
        vertex = center + center - vertex;
    }
}

void Polygon::reflex(const Line& axis) {

    for (Point& vertex: vertices) {
        double distance = axis.distanceToPoint(vertex);
        Point helpPoint(axis.a, axis.b);
        helpPoint *= (distance / (sqrt(helpPoint.x * helpPoint.x + helpPoint.y * helpPoint.y)));

        Point firstCase = vertex + helpPoint * 2.;
        Point secondCase = vertex - helpPoint * 2.;

        vertex = areEqual(axis.distanceToPoint(firstCase), distance) ? firstCase : secondCase;
    }
}

void Polygon::scale(const Point& center, const double& coefficient) {
    for (Point& vertex: vertices) {
        double newX = coefficient * (vertex.x - center.x) + center.x;
        double newY = coefficient * (vertex.y - center.y) + center.y;
        vertex = Point(newX, newY);
    }
}





/////////////////RECTANGLE CLASS/////////////////
class Rectangle : public Polygon {
public:
    Rectangle(const Point&, const Point&, double);
    Rectangle() : recCenter(Point(0.0, 0.0)), longSide(0.), shortSide(0.) {}

    Point center() const;
    std::pair<Line, Line> diagonals() const;

    double perimeter() const override;
    double area() const override;

protected:
    Point recCenter;
    double longSide;
    double shortSide;
};

Rectangle::Rectangle(const Point& first, const Point& second, double proportion) {
    double diagonalLen = getSqrtDistance(first, second);
    proportion = std::max(proportion, 1./proportion);

    recCenter = first + (second - first) * 0.5;;
    shortSide = diagonalLen / sqrt(1. + proportion * proportion);
    longSide = shortSide * proportion;

    Line sideAC(first, second);
    Point point(sideAC.a, sideAC.b);
    Point vectorAC = second - first;
    Point vectorAH = vectorAC  / (1. + proportion * proportion);
    Point H = first + vectorAH;
    point = point / sqrt(point.x * point.x + point.y * point.y);

    Point firstCase = H + point * ((getSqrtDistance(first, second) * proportion) / (proportion * proportion + 1));
    Point secondCase = H - point * ((getSqrtDistance(first, second) * proportion) / (proportion * proportion + 1));

    vertices.push_back(first);
    vertices.push_back((areEqual(getDistance(firstCase, first), shortSide) &&
                   crossProduct(point, vectorAC) > EPSILON) ? firstCase : secondCase);
    vertices.push_back(second);
    vertices.push_back(vertices[1] + (recCenter - vertices[1]) * 2);
    verticesNum = 4;
}

Point Rectangle::center() const{
    return recCenter;
}

std::pair<Line, Line> Rectangle::diagonals() const {
    return std::make_pair(Line(vertices[0], vertices[2]),
                          Line(vertices[1], vertices[3]));
}

double Rectangle::perimeter() const {
    return longSide * 2. + shortSide * 2.;
}

double Rectangle::area() const {
    return longSide * shortSide;
}




/////////////////SQUARE CLASS/////////////////
class Square : public Rectangle {
public:
    Square(const Point&, const Point&);

    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;

};

Square::Square(const Point& first, const Point& second) {
    Point shift = first - second;
    Point middle = (first + second) * 0.5;
    verticesNum = 4;
    vertices.push_back(first);
    vertices.push_back(Point(middle.x - shift.y * 0.5, middle.y + shift.x * 0.5));
    vertices.push_back(second);
    vertices.push_back(Point(middle.x + shift.y * 0.5, middle.y - shift.x * 0.5));
    longSide = getSqrtDistance(vertices[0], vertices[1]);
    shortSide = longSide;
    recCenter = middle;
}

Circle Square::inscribedCircle() const {
    return Circle(recCenter, shortSide * 0.5);
}

Circle Square::circumscribedCircle() const {
    return Circle(recCenter, 0.5 * getSqrtDistance(vertices[0], vertices[2]));
}




/////////////////TRIANGLE CLASS/////////////////
class Triangle : public Polygon {
public:
    Triangle(const Point& first, const Point& second, const Point& third) : Polygon(first, second, third) {};

    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
    Point centroid() const;
    Point orthocenter() const;
    Line EulerLine() const;
    Circle ninePointsCircle() const;

};

Circle Triangle::circumscribedCircle() const{
    double ABlen = getSqrtDistance(vertices[0], vertices[1]);
    double BClen = getSqrtDistance(vertices[1], vertices[2]);
    double CAlen = getSqrtDistance(vertices[2], vertices[0]);

    double firstNorm = vertices[0].x * vertices[0].x + vertices[0].y * vertices[0].y;
    double secondNorm = vertices[1].x * vertices[1].x + vertices[1].y * vertices[1].y;
    double thirdNorm = vertices[2].x * vertices[2].x + vertices[2].y * vertices[2].y;

    double det = 2 *(vertices[0].x * (vertices[1].y - vertices[2].y) + vertices[1].x * (vertices[2].y - vertices[0].y) +
                 vertices[2].x * (vertices[0].y - vertices[1].y));

    double x = -1. * (vertices[0].y * (secondNorm - thirdNorm) + vertices[1].y * (thirdNorm - firstNorm) +
                vertices[2].y * (firstNorm - secondNorm)) / det;
    double y = (vertices[0].x * (secondNorm - thirdNorm) + vertices[1].x * (thirdNorm - firstNorm) +
                vertices[2].x * (firstNorm - secondNorm)) / det;

    double radius = (ABlen * BClen * CAlen) / (4. * area());
    return Circle(Point(x, y), radius);
}

Circle Triangle::inscribedCircle() const{
    double ABlen = getSqrtDistance(vertices[0], vertices[1]);
    double BClen = getSqrtDistance(vertices[1], vertices[2]);
    double CAlen = getSqrtDistance(vertices[2], vertices[0]);
    Point center = (BClen * vertices[0] + CAlen * vertices[1] + ABlen * vertices[2]) / perimeter();
    return Circle(center, (2. * area()) / perimeter());
}

Point Triangle::centroid() const{
    return (vertices[0] + vertices[1] + vertices[2]) / 3;
}

Point Triangle::orthocenter() const {
    double x0 = vertices[0].x;
    double y0 = vertices[0].y;

    double x1 = vertices[1].x;
    double y1 = vertices[1].y;

    double x2 = vertices[2].x;
    double y2 = vertices[2].y;

    double firstCentPoint = y0 * ((x2 * x0 + y1 * y1) - (x0 * x1 + y2 * y2)) -
                            (x1 * x2 + y0 * y0) * (y1 - y2) +
                            ((y1 * (x0 * x1 + y2 * y2) -
                            y2 * (x2 * x0 + y1 * y1)));

    double secondCentPoint = (x0 * x0 + y1 * y2) * (x1 - x2) -
                             x0 * (x1 * x1 + y2 * y0 - (x2 * x2 + y0 * y1)) +
                             (x2 * (x1 * x1 + y2 * y0) -
                             x1 * (x2 * x2 + y0 * y1));

    double det = x0 * (y1 - y2) -
                 y0 * (x1 - x2) +
                 (x1 * y2 - x2 * y1);

    double x = firstCentPoint / det;
    double y = secondCentPoint / det;

    return Point(x, y);
}

Line Triangle::EulerLine() const{
    return {circumscribedCircle().center(), orthocenter()};
}

Circle Triangle::ninePointsCircle() const{
    return Circle((orthocenter() + circumscribedCircle().center()) * 0.5, circumscribedCircle().radius() * 0.5);
}