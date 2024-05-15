import math

import math

class Point:
    x_factor = 1.0
    y_factor = 1.0

    def __init__(self, *args, isPositive=True):
        self.x = 0.0
        self.y = 0.0
        self.isPositive = isPositive

        if len(args) == 1:
            if isinstance(args[0], Point):
                self.set(args[0].x, args[0].y, self.isPositive)
            elif isinstance(args[0], (int, float)):
                self.set(args[0], 0.0, self.isPositive)
        elif len(args) == 2:
            if isinstance(args[0], Point):
                self.set(args[0].x + args[1], args[0].y + args[1], self.isPositive)
            elif all(isinstance(arg, (int, float)) for arg in args):
                self.set(*args, self.isPositive)
        elif len(args) == 3:
            if isinstance(args[0], Point) and all(isinstance(arg, (int, float)) for arg in args[1:]):
                self.set(args[0].x + args[1], args[0].y + args[2], self.isPositive)

    def set(self, x, y, isPositive=True):
        if isPositive:
            self.x = max(x, 0)
            self.y = max(y, 0)
        else:
            self.x = x
            self.y = y

    def setX(self, x, isPositive=None):
        if isPositive is None:
            isPositive = self.isPositive
        self.set(x, self.y, isPositive)

    def setY(self, y, isPositive=None):
        if isPositive is None:
            isPositive = self.isPositive
        self.set(self.x, y, isPositive)

    def offset(self, *args):
        if len(args) == 1 and isinstance(args[0], Point):
            self.set(self.x + args[0].x, self.y + args[0].y, self.isPositive)
        elif len(args) == 2:
            self.set(self.x + args[0], self.y + args[1], self.isPositive)

    def midPoint(self, other):
        return Point((self.x + other.x) / 2, (self.y + other.y) / 2, isPositive=self.isPositive)

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y, isPositive=False)

    def __mul__(self, other):
        if isinstance(other, Point):
            return self.x * other.x + self.y * other.y
        elif isinstance(other, (int, float)):
            return Point(self.x * other, self.y * other, isPositive=False)

    def __rmul__(self, other):
        return self.__mul__(other)

    def expand(self):
        self.setX(self.x * Point.x_factor)
        self.setY(self.y * Point.y_factor)

    @classmethod
    def set_x_factor(cls, x_factor):
        cls.x_factor = x_factor

    @classmethod
    def set_y_factor(cls, y_factor):
        cls.y_factor = y_factor

    @classmethod
    def set_x_y_factors(cls, x_factor, y_factor):
        cls.x_factor = x_factor
        cls.y_factor = y_factor

    def dist_to_line(self, line):
        if line.slope == float('inf'):
            return abs(self.x - (-line.intercept))
        return abs(line.slope * self.x - self.y + line.intercept) / math.sqrt(line.slope ** 2 + 1)

    def __str__(self):
        return f"Point(x={self.x}, y={self.y})"

    def display(self):
        print(self.__str__())

    def __eq__(self, other):
        if isinstance(other, Point):
            return self.x == other.x and self.y == other.y
        return False





class LineSegment:
    def __init__(self, *args):
        # Constructor for two Points
        if len(args) == 2 and all(isinstance(arg, Point) for arg in args):
            self.start = args[0]
            self.end = args[1]
        # Constructor for copy of another LineSegment
        elif len(args) == 1 and isinstance(args[0], LineSegment):
            self.start = Point(args[0].start)
            self.end = Point(args[0].end)
        # Constructor for middle point, length, and slope
        elif len(args) == 3 and isinstance(args[0], Point) and all(isinstance(arg, (int, float)) for arg in args[1:]):
            middle, length, slope = args
            dx = (length / 2) * math.cos(math.atan(slope))
            dy = (length / 2) * math.sin(math.atan(slope))
            self.start = Point(middle.x - dx, middle.y - dy, isPositive=True)
            self.end = Point(middle.x + dx, middle.y + dy, isPositive=True)

    def length(self):
        return math.sqrt((self.end.x - self.start.x)**2 + (self.end.y - self.start.y)**2)

    def midPoint(self):
        mid_x = (self.start.x + self.end.x) / 2
        mid_y = (self.start.y + self.end.y) / 2
        return Point(mid_x, mid_y)

    def intersection(self, other):
        # Convert to line objects to find the intersection
        line1 = Line(self.start, self.end)
        line2 = Line(other.start, other.end)
        intersection_point = line1.intersectLine(line2)

        if intersection_point and self.isOnLineSegment(intersection_point) and other.isOnLineSegment(intersection_point):
            return intersection_point
        return None

    def isOnLineSegment(self, point):
        # Check if the point is within the bounding box defined by the segment's endpoints
        within_x = min(self.start.x, self.end.x) <= point.x <= max(self.start.x, self.end.x)
        within_y = min(self.start.y, self.end.y) <= point.y <= max(self.start.y, self.end.y)
        return within_x and within_y

    def __str__(self):
        return f"LineSegment(start={self.start}, end={self.end})"

    def display(self):
        print(self.__str__())

    def __eq__(self, other):
        if isinstance(other, LineSegment):
            return (self.start == other.start and self.end == other.end) or (self.start == other.end and self.end == other.start)
        return False




class Line:
    def __init__(self, *args):
        # Different constructors based on input arguments
        if len(args) == 2 and all(isinstance(arg, (int, float)) for arg in args):
            self.slope = args[0]
            self.intercept = args[1]
        elif len(args) == 3 and all(isinstance(arg, (int, float)) for arg in args):
            a, b, c = args
            if b != 0:
                self.slope = -a / b
                self.intercept = c / b
            else:
                self.slope = float('inf')  # Infinite slope represents a vertical line
                self.intercept = None
        elif len(args) == 2 and all(isinstance(arg, Point) for arg in args):
            p1, p2 = args
            if p1.x != p2.x:
                self.slope = (p2.y - p1.y) / (p2.x - p1.x)
                self.intercept = p1.y - self.slope * p1.x
            else:
                self.slope = float('inf')
                self.intercept = None
        elif len(args) == 1 and isinstance(args[0], LineSegment):
            p1, p2 = args[0].start, args[0].end
            if p1.x != p2.x:
                self.slope = (p2.y - p1.y) / (p2.x - p1.x)
                self.intercept = p1.y - self.slope * p1.x
            else:
                self.slope = float('inf')
                self.intercept = None
        elif len(args) == 2 and isinstance(args[0], Point) and isinstance(args[1], (int, float)):
            p, slope = args
            self.slope = slope
            self.intercept = p.y - slope * p.x

    def intersectLine(self, other):
        if self.slope == other.slope:
            if self.intercept == other.intercept:
                return None  # Lines are the same
            else:
                return None  # Parallel lines do not intersect
        elif self.slope == float('inf'):  # self is vertical
            x = -self.intercept
            y = other.slope * x + other.intercept
            return Point(x, y, isPositive=False)
        elif other.slope == float('inf'):  # other is vertical
            x = -other.intercept
            y = self.slope * x + self.intercept
            return Point(x, y, isPositive=False)
        else:
            x = (other.intercept - self.intercept) / (self.slope - other.slope)
            y = self.slope * x + self.intercept
            return Point(x, y, isPositive=False)

    def distance(self, point):
        if self.slope == float('inf'):
            return abs(point.x + self.intercept)
        else:
            return abs(self.slope * point.x - point.y + self.intercept) / ((self.slope**2 + 1)**0.5)

    def perpendicular(self, point):
        if self.slope == 0:
            # Vertical line
            return Line(point, float('inf'))
        elif self.slope == float('inf'):
            # Horizontal line
            return Line(point, 0)
        else:
            # Perpendicular slope is the negative reciprocal
            new_slope = -1 / self.slope
            return Line(point, new_slope)

    def __str__(self):
        if self.slope == float('inf'):
            return f"x = {-self.intercept}"
        else:
            return f"y = {self.slope}x + {self.intercept}"

    def display(self):
        print(self.__str__())

    def __eq__(self, other):
        return (self.slope == other.slope) and (self.intercept == other.intercept)




class Triangle:
    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C

    def area(self):
        # Using the shoelace formula for area of a triangle on a coordinate plane
        return abs((self.A.x * (self.B.y - self.C.y) + self.B.x * (self.C.y - self.A.y) + self.C.x * (self.A.y - self.B.y)) / 2)

    def perimeter(self):
        # Perimeter is the sum of lengths of all sides
        sideAB = math.sqrt((self.B.x - self.A.x)**2 + (self.B.y - self.A.y)**2)
        sideBC = math.sqrt((self.C.x - self.B.x)**2 + (self.C.y - self.B.y)**2)
        sideCA = math.sqrt((self.A.x - self.C.x)**2 + (self.A.y - self.C.y)**2)
        return sideAB + sideBC + sideCA

    def angles(self):
        # Calculate angles using the cosine rule
        def angle(a, b, c):
            # Calculate angle opposite side c using cosine rule
            cos_theta = (a**2 + b**2 - c**2) / (2 * a * b)
            return math.acos(cos_theta) * (180 / math.pi)  # Converting radians to degrees
        
        sideAB = math.sqrt((self.B.x - self.A.x)**2 + (self.B.y - self.A.y)**2)
        sideBC = math.sqrt((self.C.x - self.B.x)**2 + (self.C.y - self.B.y)**2)
        sideCA = math.sqrt((self.A.x - self.C.x)**2 + (self.A.y - self.C.y)**2)

        angleA = angle(sideBC, sideCA, sideAB)
        angleB = angle(sideCA, sideAB, sideBC)
        angleC = angle(sideAB, sideBC, sideCA)

        return (angleA, angleB, angleC)

    def circumscribedCircle(self):
        # Circumscribed circle (circumcircle) center and radius
        # Using perpendicular bisectors of sides AB and BC
        def midpoint(p1, p2):
            return Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2)

        def perpendicular_bisector(p1, p2):
            m = midpoint(p1, p2)
            if p2.x != p1.x:
                slope = -(p2.x - p1.x) / (p2.y - p1.y)
                intercept = m.y - slope * m.x
                return (slope, intercept)
            else:
                return None  # vertical line

        def intersection(line1, line2):
            if line1 and line2:
                # Normal lines
                slope1, intercept1 = line1
                slope2, intercept2 = line2
                x = (intercept2 - intercept1) / (slope1 - slope2)
                y = slope1 * x + intercept1
                return Point(x, y)
            elif line1 is None:
                # Line1 vertical
                _, intercept2 = line2
                x = midpoint(self.A, self.B).x
                y = line2[0] * x + intercept2
                return Point(x, y)
            elif line2 is None:
                # Line2 vertical
                _, intercept1 = line1
                x = midpoint(self.B, self.C).x
                y = line1[0] * x + intercept1
                return Point(x, y)

        bisectorAB = perpendicular_bisector(self.A, self.B)
        bisectorBC = perpendicular_bisector(self.B, self.C)
        center = intersection(bisectorAB, bisectorBC)
        radius = math.sqrt((self.A.x - center.x)**2 + (self.A.y - center.y)**2)

        return (center, radius)

    def __str__(self):
        return f"Triangle(A={self.A}, B={self.B}, C={self.C})"

    def display(self):
        print(self.__str__())

    def __eq__(self, other):
        # Checking if two triangles are the same, irrespective of the order of points
        return sorted([self.A, self.B, self.C]) == sorted([other.A, other.B, other.C])

# Assuming Point, Line, LineSegment, and Triangle classes are already defined as provided above

def main():
    # Create Points
    p1 = Point(3, 4)
    p2 = Point(1, 2)
    p3 = Point(5, 1)
    p4 = Point(-1, -1, isPositive=False)  # Demonstrating a point allowed in any quadrant
    p5 = Point(p1)  # Copy of p1

    # Display Points
    p1.display()
    p2.display()
    p3.display()
    p4.display()
    p5.display()

    # Create Line Segments
    seg1 = LineSegment(p1, p2)
    seg2 = LineSegment(p2, p3)
    seg3 = LineSegment(p3, p1)

    # Display Line Segments and their lengths
    seg1.display()
    seg2.display()
    seg3.display()
    print(f"Length of seg1: {seg1.length()}")
    print(f"Length of seg2: {seg2.length()}")
    print(f"Length of seg3: {seg3.length()}")

    # Create Lines from Points
    line1 = Line(p1, p2)
    line2 = Line(p3, p4)

    # Display Lines
    line1.display()
    line2.display()

    # Intersect Lines
    intersection_point = line1.intersectLine(line2)
    if intersection_point:
        print("Intersection Point:")
        intersection_point.display()
    else:
        print("Lines do not intersect.")

    # Create a Triangle
    triangle = Triangle(p1, p2, p3)
    triangle.display()

    # Display Triangle Properties
    print(f"Area of Triangle: {triangle.area()}")
    print(f"Perimeter of Triangle: {triangle.perimeter()}")
    angles = triangle.angles()
    print(f"Angles of Triangle: {angles[0]:.2f}, {angles[1]:.2f}, {angles[2]:.2f}")
    center, radius = triangle.circumscribedCircle()
    print(f"Circumscribed Circle Center: ({center.x}, {center.y}), Radius: {radius:.2f}")

    # Modify a point and see the effect on Line Segment
    print("Modifying p1...")
    p1.setX(10)
    p1.setY(10)
    seg1.display()  # Showing how the line segment changes when one of its points is changed

if __name__ == "__main__":
    main()

