using Autodesk.AutoCAD.ApplicationServices;
using Autodesk.AutoCAD.DatabaseServices;
using Autodesk.AutoCAD.EditorInput;
using Autodesk.AutoCAD.Geometry;
using Autodesk.AutoCAD.Runtime;

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MgdAcApplication = Autodesk.AutoCAD.ApplicationServices.Application;
using MgdAcDocument = Autodesk.AutoCAD.ApplicationServices.Document;
using AcWindowsNS = Autodesk.AutoCAD.Windows;

namespace PatternWalker
{
    //KNOWN BUGS
    //POLYLINES have to be closed for code to find the option distance paths
    //program seems to freeze on certain parts
    public class MutateJig : DrawJig
    {
        #region fields
        //declare internal state
        public double m_baseSlope, m_WalkSlope, m_rotAngle;
        public Polyline m_pLineAnchor, m_pLineWalk;
        public Point3d m_anchorPoint, m_location, m_walkPoint, m_walkBase;
        private bool m_plineAnchorCW, m_plineWalkCW, m_lsegAnchorCW;
        #endregion

        //Constructor sets the state
        //*if a loop adds all selection of ents to a list, can scale to jig complex parts
        public MutateJig(Polyline pLineAnchor, Polyline pLineWalk, Point3d anchorPoint, Point3d walkBasePoint)
        {
            m_anchorPoint = anchorPoint.TransformBy(UCS);
            m_pLineAnchor = pLineAnchor;
            m_pLineWalk = pLineWalk;
            m_walkBase = walkBasePoint.TransformBy(UCS);
        }

        #region Properties
        public void TransformPolyline()
        {
            Matrix3d mat = Transformation;
            m_pLineWalk.TransformBy(mat);
        }

        public double Angle
        {
            get { return m_rotAngle; }
            set { m_rotAngle = value; }
        }

        public double anchorSlope
        {
            get { return m_baseSlope; }
            set { m_baseSlope = value; }
        }

        public double walkSlope
        {
            get { return m_WalkSlope; }
            set { m_WalkSlope = value; }
        }

        public Polyline Anchor
        {
            get { return m_pLineAnchor; }
            set { m_pLineAnchor = value; }
        }

        public Polyline Walker
        {
            get { return m_pLineWalk; }
            set { m_pLineWalk = value; }
        }

        public Point3d anchBasePoint
        {
            get { return m_anchorPoint; }
            set { m_anchorPoint = value; }
        }

        public Point3d ConnectedPoint
        {
            get { return m_location; }
            set { m_location = value; }
        }

        public Point3d walkpoint
        {
            get { return m_walkPoint; }
            set { m_walkPoint = value; }
        }

        public Editor AcadEd
        {
            get
            { return MgdAcApplication.DocumentManager.MdiActiveDocument.Editor; }
        }

        public Matrix3d UCS
        {
            get { return AcadEd.CurrentUserCoordinateSystem; }
        }

        public Matrix3d Transformation
        {
            get
            {
                return
                Matrix3d.Rotation(m_rotAngle, Vector3d.ZAxis.TransformBy(UCS), m_location).
                PostMultiplyBy(Matrix3d.Displacement(m_walkPoint.GetVectorTo(m_location)));
            }
        }
        #endregion

        #region Overrides
        protected override SamplerStatus Sampler(JigPrompts prompts)
        {
            //aquire a point on the anchored polyline

            JigPromptPointOptions jpo = new JigPromptPointOptions("\nPoint to connect patterns: ");
            jpo.UserInputControls = UserInputControls.GovernedByOrthoMode | UserInputControls.GovernedByUCSDetect;
            PromptPointResult ppr = prompts.AcquirePoint(jpo);

            if (ppr.Status == PromptStatus.Cancel && ppr.Status == PromptStatus.Error)
                return SamplerStatus.Cancel;

            if (ppr.Status == PromptStatus.OK)
            {
                //if no change, don't do anything *bc of flicker
                if (m_anchorPoint == ppr.Value)
                { return SamplerStatus.NoChange; }
                else
                {
                    //set point to be the closest to the pline
                    m_location = ppr.Value;
                    m_location = m_pLineAnchor.GetClosestPointTo(m_location, true);

                    double distance = findDistance(m_anchorPoint, m_pLineAnchor);
                    m_walkPoint = findWalkPoint(m_walkBase, m_pLineWalk, !m_lsegAnchorCW, distance);

                    //convert slopes to angles (rad)
                    double anchorAng = slopeToAngle(m_baseSlope);
                    double walkAng = slopeToAngle(m_WalkSlope);

                    m_rotAngle = anchorAng - walkAng + Math.PI;
                    //m_rotAngle = Math.Atan((m_baseSlope - m_WalkSlope) / (1 + m_baseSlope * m_WalkSlope));

                    m_plineAnchorCW = clockwise(m_pLineAnchor);
                    m_plineWalkCW = clockwise(m_pLineWalk);

                    return SamplerStatus.OK;
                }
            }
            return SamplerStatus.Cancel;
        }

        protected override bool WorldDraw(Autodesk.AutoCAD.GraphicsInterface.WorldDraw draw)
        {
            Matrix3d trans = Transformation;
            Autodesk.AutoCAD.GraphicsInterface.WorldGeometry geo = draw.Geometry;

            if (geo != null)
            {
                geo.PushModelTransform(trans);
                geo.Draw(m_pLineWalk);
                geo.PopModelTransform();
            }
            return true;
        }

        #endregion

        //convert the slope to an angle from the x-axis moving counterclockwise
        public double slopeToAngle(double m)
        {
            double angle;                 
            //if the slope is in and infinite state (devided by 0) then we know it is vertical line
            if (m == double.NegativeInfinity || m == double.PositiveInfinity)
            { 
                angle = Math.PI / 2;
                return angle;
            }

            angle = Math.Atan(m); 
            return angle;
        }

        public Point3d findWalkPoint(Point3d basePoint, Polyline pLine,
        bool CW,
        double distance)
        {
            bool plineCw = clockwise(pLine);
            int direction = 0;
            //the Clockwise bool being passed in is for the direction we want our line to go
            //if the second poly line is clockwise then it directly corallates to the line seg
            //if the second poly line is CCW then it inversely corallates to the line segment
            if (plineCw == true)
            {
                if (CW == true)
                { direction = 1; }
                else
                { direction = -1; }
            }
            else
            {
                if (CW == true)
                { direction = -1; }
                else
                { direction = 1; }
            }

            //find the part the basePoint is on
            int startPart = findPart(basePoint, pLine);
            //use direction to decide if ascending or desceding
            //get first segment dist from base point to start or end *based on direction
            double runningDist = segmentDist(basePoint, pLine, startPart, direction);
            SegmentType segType = new SegmentType();
            int index = startPart;
            Curve3d seg = null;

            //if base point isnt on a vertex, need to get the piece of the segment first, then add it accordingly

            //add the following segments until distance is greater than distance being sought
            while (runningDist < distance)
            {
                index = index + direction;
                //make sure the index doesnt become an out of bounds vertex
                if (index == -1)
                { index = pLine.NumberOfVertices - 1; }
                else if (index == pLine.NumberOfVertices)
                { index = 0; }

                //get the segment of each piece and add it on
                segType = pLine.GetSegmentType(index);
                if (segType == SegmentType.Line)
                {
                    LineSegment3d lseg = pLine.GetLineSegmentAt(index);
                    runningDist = runningDist + lseg.Length;
                }
                else
                {
                    CircularArc3d aSeg = pLine.GetArcSegmentAt(index);
                    double r = aSeg.Radius;
                    double liDist = aSeg.StartPoint.DistanceTo(aSeg.EndPoint);
                    double dist = (.5 * liDist) / r;
                    dist = r * 2 * Math.Asin(dist);

                    runningDist = runningDist + dist;
                }
            }

            //get the segment used last
            segType = pLine.GetSegmentType(index);
            if (segType == SegmentType.Line)
            { seg = pLine.GetLineSegmentAt(index); }
            else
            { seg = pLine.GetArcSegmentAt(index); }

            //determine the last point in segment that both the distance and running distance would have
            Point3d sharedPoint = new Point3d();
            Point3d lastPoint = new Point3d();
            //depends on direction
            if (direction > 0)
            {
                sharedPoint = seg.StartPoint;
                lastPoint = seg.EndPoint;
            }
            else if (direction < 0)
            { 
                sharedPoint = seg.EndPoint;
                lastPoint = seg.StartPoint;
            }

            //find the point on the last segment using the excess distance between the two
            double remainder = Math.Abs(runningDist - distance);
            Point3d finalPoint = new Point3d();
            //if remainder is 0 then it was on an endpoint we ended
            if (remainder == 0)
            {
                finalPoint = lastPoint;
                m_WalkSlope = slopeAtPoint(finalPoint, m_pLineWalk, index);
                return lastPoint;
            }
            //otherwise need to get the last segment distance and find the point as it would fall               
            else
            {
                //back up from the the endpoint of the total distance of the last segment the distance of the remainder
                if (segType == SegmentType.Line)
                {
                    LineSegment3d lDat = pLine.GetLineSegmentAt(index);
                    List<Point3d> points = lineCircIntersect(lDat, remainder, lastPoint);
                    //verify that the points are on the segment
                    points = isPointOnCurve(points, pLine);
                    //determine which point is closer to the previous point (shared Point)
                    finalPoint = findClosestPoint(points, sharedPoint);
                    
                    m_WalkSlope = slopeAtPoint(finalPoint, m_pLineWalk, index);
                    return finalPoint;
                }
                else
                {
                    CircularArc3d aDat = pLine.GetArcSegmentAt(index);
                    List<Point3d> points = notchOnArc(aDat, remainder, lastPoint);
                    //verify that the points are on the segment
                    points = isPointOnCurve(points, pLine);
                    //determine which point is closer to the previous point
                    finalPoint = findClosestPoint(points, sharedPoint);
                    m_WalkSlope = slopeAtPoint(finalPoint, m_pLineWalk, index);
                    return finalPoint;
                }
            }

            //this point will be used to calculate slopes *for rotation and rotation point
        }


        //calculates intersections between a line and circle
        public List<Point3d> lineCircIntersect(LineSegment3d lDat, double r, Point3d center)
        {
            //find known constants
            double m = getSlopeLine(lDat);
            double p = center.X;
            double q = center.Y;
            double c = findYint(lDat.StartPoint, m);

            List<Point3d> points = new List<Point3d>();
            Point3d first;
            Point3d second;

            //if slope is 0 then it is horizontal
            if (m == 0)
            {
                first = new Point3d(center.X + r, center.Y, center.Z);
                second = new Point3d(center.X - r, center.Y, center.Z);
            }
            //if slope is an infinity value it is vertical
            else if (m == double.PositiveInfinity || m == double.NegativeInfinity)
            {
                first = new Point3d(center.X, center.Y + r, center.Z);
                second = new Point3d(center.X, center.Y - r, center.Z);
            }
            else
            {
                //using line equation  y=mx+c
                //using circle equation (x-p)^2 + (y-q)^2 =r^2   (p,h) = center
                //----------------
                // (x-p)^2 + (mx+c-q)^2 = r^2
                // (m^2 +1)x^2 = 2(mc - mq-p)x + (q^2 - r^2 + p^2 - 2cq + c^2) = 0
                //     A              B                         C
                //-----
                //quadratic formula  x= (-B +/- (B2 - 4AC)^.5) / 2A
                // plug back into line
                //   y= m ((-B +/- (B2 - 4AC)^.5) / 2A) + c
                // plug A,B,C back into the formula

                double quadNumerator = -2 * (m * c - m * q - p) + Math.Sqrt(Math.Pow((2 * (m * c - m * q - p)), 2) - 4 * (m * m + 1) * (q * q - r * r + p * p - 2 * c * q + c * c));
                double quadNumNeg = -2 * (m * c - m * q - p) - Math.Sqrt(Math.Pow((2 * (m * c - m * q - p)), 2) - 4 * (m * m + 1) * (q * q - r * r + p * p - 2 * c * q + c * c));
                double quadDenom = 2 * (m * m + 1);
                double QuadForm = quadNumerator / quadDenom;
                double QuadFormNeg = quadNumNeg / quadDenom;
                double y = m * QuadForm + c;
                double yNeg = m * QuadFormNeg + c;

                //use the +/- results to solve for X
                double x = (y - c) / m;
                double xNeg = (yNeg - c) / m;

                first = new Point3d(x, y, 0);
                second = new Point3d(xNeg, yNeg, 0);
            }

            //return the two points            
            points.Add(first);
            points.Add(second);

            return points;
        }

        public List<Point3d> notchOnArc(CircularArc3d aDat, double remainder, Point3d sharedPoint)
        {
            List<Point3d> points = new List<Point3d>();

            Point3d insertPoint1 = new Point3d();
            Point3d insertPoint2 = new Point3d();

            Point3d center = aDat.Center;
            double alpha = remainder / aDat.Radius;
            double beta = Math.Acos(((sharedPoint.X - aDat.Center.X) / aDat.Radius));

            //grabbing the angle from CAD ignored which point we were measuring from
            //if (aDat.StartAngle > aDat.EndAngle)
            //{ beta = aDat.StartAngle; }
            //else
            //{ beta = aDat.EndAngle; }

            double theta = beta - alpha;
            double xDim = aDat.Radius * Math.Round(Math.Cos(theta), 9);
            double yDim = aDat.Radius * Math.Round(Math.Sin(theta), 9);

            double xDimAlt = aDat.Radius * Math.Round(Math.Cos(beta + alpha), 9);
            double yDimAlt = aDat.Radius * Math.Round(Math.Sin(beta + alpha), 9);

            insertPoint1 = new Point3d(center.X + xDim, center.Y + yDim, 0);
            insertPoint2 = new Point3d(center.X + xDimAlt, center.Y + yDimAlt, 0);

            points.Add(insertPoint1);
            points.Add(insertPoint2);

            return points;
        }

        //get the y intersect
        private double findYint(Point3d point, double m)
        {
            double yIntercept = point.Y - m * point.X;
            return yIntercept;
        }

        //flaw to this is that its possible that both points exist on the polyline,
        //needs to look specifically at just the part segment, not whole polyline
        //might simply have it loop through the poly line and make DBObjects that
        //can be used to compare to the endpoint,startpoint of the seg data
        //if one matches then it gets passed and used for a Curve
        public List<Point3d> isPointOnCurve(List<Point3d> points, Polyline pLine)
        {
            List<Point3d> returnPoints = new List<Point3d>();
            foreach (Point3d p in points)
            {
                try
                {
                    Point3d pt = pLine.GetClosestPointTo(p, false);
                    double d = (pt - p).Length;
                    if (d < 0.001)
                    { returnPoints.Add(p); }
                }
                catch { }
            }

            return returnPoints;
        }

        //hand over the two intersection points, calculate distance from the previous endPoint
        //whichever is the closest point gets returned
        public Point3d findClosestPoint(List<Point3d> points, Point3d endPoint)
        {
            double distance = 0;
            Point3d returnMe = new Point3d();
            foreach (Point3d p in points)
            {
                double d = Math.Sqrt((Math.Pow((p.X - endPoint.X), 2) + Math.Pow((p.Y - endPoint.Y), 2)));
                if (d < distance || distance == 0)
                {
                    returnMe = p;
                    distance = d;
                }
            }
            return returnMe;
        }

        //vvvvvvvvvvvvvvvv need to pass it the received point from drag
        public double findDistance(Point3d anchor, Polyline pLine)
        {
            double distance = 0;
            //first find the part of the polyline that this anchor point is on
            int startPart = findPart(anchor, pLine);
            // -1 is returned if point isn't found
            if (startPart == -1)
            { return 0; }

            //need to get a point for the end to pass for the last part
            //point needs to be point nearest on the poly line from the curser
            Point3d nextPoint = m_location;

            int lastPart = findPart(nextPoint, pLine);
            m_baseSlope = slopeAtPoint(nextPoint, m_pLineAnchor, lastPart);

            //if the points are on the same part we need to check distance
            #region same part
            if (lastPart == startPart)
            {
                Curve3d seg = null;
                SegmentType segType = pLine.GetSegmentType(lastPart);
                if (segType == SegmentType.Line)
                {
                    seg = pLine.GetLineSegmentAt(lastPart);
                    LineSegment3d lDat = pLine.GetLineSegmentAt(lastPart);
                    //dist form
                    distance = anchor.DistanceTo(nextPoint);
                }
                else
                {
                    seg = pLine.GetArcSegmentAt(lastPart);
                    CircularArc3d aDat = pLine.GetArcSegmentAt(lastPart);
                    double r = aDat.Radius;
                    double liDist = anchor.DistanceTo(nextPoint);
                    distance = (.5 * liDist) / r;
                    distance = r * 2 * Math.Asin(distance);
                }


                
                //descern which points are closet to which end
                Point3d vertexPnt = pLine.GetPoint3dAt(lastPart);
                double anchorToStart = anchor.DistanceTo(seg.StartPoint);
                double nextToStart = nextPoint.DistanceTo(seg.StartPoint);
                if (vertexPnt == seg.StartPoint)
                {
                    if (anchorToStart < nextToStart)
                    {
                        //then anchor point is closer to start point which is the vertex we're after
                        //if we already know its a a Clockwise motion then this means the
                        //line we just made was a Clockwise path
                        if (m_plineAnchorCW == true)
                        { m_lsegAnchorCW = true; }
                        else
                        { m_lsegAnchorCW = false; }
                    }
                    else
                    {
                        //then the next point is closer to the start point
                        //if the polyline is clockwise then we are moving  counterclockwise with this line
                        if (m_plineAnchorCW == true)
                        { m_lsegAnchorCW = false; }
                        else
                        { m_lsegAnchorCW = true; }
                    }
                }

                //by designating the  point to reference the vertex # we can decide
                //if the movement was ascending or descending for CW or CCW

            #endregion
            }
            //if the points are on seperate parts
            else
            {
                //use knowledge of CW and CCW to know ahead of time which direction we're going
                double dist1 = longDistance(anchor, nextPoint, pLine, startPart, lastPart, 1);
                double dist2 = longDistance(anchor, nextPoint, pLine, startPart, lastPart, -1);
                //we'll use the shorter of the two distances <--??

                if (dist1 > dist2)
                {
                    distance = dist2;
                    if (m_plineAnchorCW == true)
                    {
                        m_lsegAnchorCW = false;
                    }
                    else
                        m_lsegAnchorCW = true;
                }
                else
                {
                    distance = dist1;
                    if (m_plineAnchorCW == true)
                    {
                        m_lsegAnchorCW = true;
                    }
                    else
                        m_lsegAnchorCW = false;

                }
                //ascending and descending should indicate ClockWise or Counter-ClockWise for the two distances
            }

            //default return, similar to null
            return distance;
        }

        //pass it the polyline, the basepoint, the targetpoint, and direction(ascending or descending)
        //should return the distance between those points along that path
        public double longDistance(Point3d basePoint, Point3d pickedPoint, Polyline pLine, int startPart, int lastPart, int direction)
        {
            double distance = 0;
            SegmentType segType = new SegmentType();

            //if direction is positive then we are moving down the vertecies, and need the endpoint
            //if the direction is negative we are moving backward along the verticies and need the startpoint
            //add in first segment
            distance = segmentDist(basePoint, pLine, startPart, direction);
            //add in last segment
            distance = distance + segmentDist(pickedPoint, pLine, lastPart, direction * -1);

            startPart = startPart + direction;

            //loop through any segments which are totally used and do not have the last point
            for (int i = startPart; i != lastPart; i = i + direction)
            {
                //make sure index doenst become a vertex which doesnt exist
                if (i == -1)
                {
                    i = pLine.NumberOfVertices - 1;
                    if (lastPart == pLine.NumberOfVertices - 1)
                    { return distance; }
                }
                else if (i == pLine.NumberOfVertices)
                {
                    i = 0;
                    if (lastPart == 0)
                    { return distance; }
                }

                segType = pLine.GetSegmentType(i);
                if (segType == SegmentType.Line)
                {
                    LineSegment3d lseg = pLine.GetLineSegmentAt(i);
                    distance = distance + lseg.Length;
                }
                else
                {
                    CircularArc3d aSeg = pLine.GetArcSegmentAt(i);
                    double r = aSeg.Radius;
                    double liDist = aSeg.StartPoint.DistanceTo(aSeg.EndPoint);
                    double dist = (.5 * liDist) / r;
                    dist = r * 2 * Math.Asin(dist);

                    distance = distance + dist;
                }
            }

            //return that distance
            return distance;
        }

        public double segmentDist(Point3d p, Polyline pLine, int partIndex, int direction)
        {
            Curve3d seg = null;
            double distance = 0;

            //if direction is positive then we are moving down the vertecies, and need the endpoint
            //if the direction is negative we are moving backward along the verticies and need the startpoint
            SegmentType segType = pLine.GetSegmentType(partIndex);
            if (segType == SegmentType.Line)
            {
                seg = pLine.GetLineSegmentAt(partIndex);
                if (direction > 0)
                    distance = seg.EndPoint.DistanceTo(p);
                else
                    distance = seg.StartPoint.DistanceTo(p);
            }
            else
            {
                seg = pLine.GetArcSegmentAt(partIndex);

                Point3d nextPoint = new Point3d();
                if (direction > 0)
                { nextPoint = seg.EndPoint; }
                else
                { nextPoint = seg.StartPoint; }

                CircularArc3d aDat = pLine.GetArcSegmentAt(partIndex);
                double r = aDat.Radius;
                double liDist = p.DistanceTo(nextPoint);
                if(liDist == 0)
                {
                    return 0;
                }
                distance = (.5 * liDist) / r;
                distance = r * 2 * Math.Asin(distance);
            }

            return distance;
        }

        public int findPart(Point3d start, Polyline pLine)
        {
            int index = -1;
            Curve3d seg = null;
            //loop through the poly line segment by segment
            for (int i = 0; i < pLine.NumberOfVertices; i++)
            {
                SegmentType segType = pLine.GetSegmentType(i);
                if (segType == SegmentType.Line)
                { seg = pLine.GetLineSegmentAt(i); }
                else
                { seg = pLine.GetArcSegmentAt(i); }

                if (seg != null)
                {
                    if (seg.IsOn(start) == true)
                        return i;
                }
            }

            return index;
        }

        public bool clockwise(Polyline pline)
        {
            bool CW = true;
            //for each vertex in the poly line
            //calculate a cross product to decide if it is pos or neg
            //if more vertices in a polygon are pos it is CCW
            //if more vertices in a polygon are neg it is CW
            int pos = 0;
            int neg = 0;

            for (int i = 0; i < pline.NumberOfVertices; i++)
            {
                double answer;
                //change the algorythm a touch for the first and last vertex
                if (i == 0)
                {
                    //pass the previous vertex, current vertex, and next vertex
                    answer = crossProduct(pline.GetPoint3dAt(pline.NumberOfVertices - 1), pline.GetPoint3dAt(i), pline.GetPoint3dAt(i + 1));
                }
                else if (i == (pline.NumberOfVertices - 1))
                { answer = crossProduct(pline.GetPoint3dAt(i - 1), pline.GetPoint3dAt(i), pline.GetPoint3dAt(0)); }
                else
                { answer = crossProduct(pline.GetPoint3dAt(i - 1), pline.GetPoint3dAt(i), pline.GetPoint3dAt(i + 1)); }

                if (answer > 0)
                { pos++; }
                else if (answer < 0)
                { neg++; }
            }

            if (pos > neg)
            { CW = false; }
            else
            { CW = true; }

            return CW;
        }

        //pass the points of the vertex we're examining and its preceding, proceding points
        //cross multiply to find the area and if protrudes pos or neg from the plane
        public double crossProduct(Point3d v1, Point3d v2, Point3d v3)
        {
            double answer = ((v2.X - v1.X) * (v3.Y - v2.Y)) - ((v2.Y - v1.Y) * (v3.X - v1.X));
            return answer;
        }

        //pass a segment, based on the slope at that point on that seg, it will use the appropriate slope algorythm
        private double slopeAtPoint(Point3d point, Polyline pLine, int index)
        {
            double m = new double();
            SegmentType segType = pLine.GetSegmentType(index);
            if (segType == SegmentType.Line)
            {
                LineSegment3d lDat = pLine.GetLineSegmentAt(index);
                m = getSlopeLine(lDat);
            }
            else
            {
                CircularArc3d aDat = pLine.GetArcSegmentAt(index);
                m = getSlopeArc(aDat, point);
            }
            return m;
        }

        private double getSlopeArc(CircularArc3d aDat, Point3d tpoint)
        {
            double slope = 0;
            double x1, x2, y1, y2;
            x2 = Math.Round(aDat.Center.X, 7);
            x1 = Math.Round(tpoint.X);
            y2 = Math.Round(aDat.Center.Y, 7);
            y1 = Math.Round(tpoint.Y, 7);

            //get slope of the radius as it touches the endpoint
            //double m = (aDat.Center.Y - tpoint.Y) / (aDat.Center.X - tpoint.X);
            double m = (y2 - y1) / (x2 - x1);

            //don't know if it has to account for undefined numbers here, code seems to account for it later
            if(m == double.NegativeInfinity | m == double.PositiveInfinity)
            { return Math.PI;}
            else if(m == 0)
            { return Math.PI/2;}
            else
                slope = -1 * (1 / m);

            //be sure to check for infinities and 0 slopes
            return slope;
        }

        private double getSlopeLine(LineSegment3d lDat)
        {
            double x1, x2, y1, y2;
            x2 = Math.Round(lDat.EndPoint.X, 7);
            x1 = Math.Round(lDat.StartPoint.X, 7);
            y2 = Math.Round(lDat.EndPoint.Y, 7);
            y1 = Math.Round(lDat.StartPoint.Y, 7);
            //double m = (Math.Round(lDat.EndPoint.Y, 9) - Math.Round(lDat.StartPoint.Y,9)) / (Math.Round(lDat.EndPoint.X,9) - Math.Round(lDat.StartPoint.Y,9));
            //need to be sure to handle infinties and 0 slopes
            double m = (y2 - y1) / (x2 - x1);
            return m;
        }
    }

    //*****TO DO
    //need to account for vertical and horizontal slopes, to be sure that neg and pos infinity dont break code
    //need to fix isPointOnCurve to look at a particular segment instead of whole polyline
    //need to solve the draw to keep all data of the jig drag, somehow it is making the changes to the original
    //  position and rotation
    //need to look through code to refine process to require much less processing, to smoothe out the flickering
    public class Class1
    { 
    [CommandMethod("PatternWalk")]
    public void start()
    {
        Document doc = Application.DocumentManager.MdiActiveDocument;
        Editor ed = doc.Editor;
        Database db = doc.Database;

        using (Transaction tr = db.TransactionManager.StartTransaction())
        {
            BlockTable bt = tr.GetObject(db.BlockTableId, OpenMode.ForRead) as BlockTable;
            //request the anchor polyline          
            string msg = "\nSelect anchor pattern (only 1 PolyLine)";
            Polyline anchorPline = getPolyLine(tr, ed, msg);
            if (anchorPline == null)
            { return; }

            PromptPointResult pntRes;
            PromptPointOptions ppo = new PromptPointOptions("");
            ppo.Message = "\nSelect start point from anchored pattern.";

            //request start point
            Point3d anchorPoint = new Point3d();
            pntRes = ed.GetPoint(ppo);
            if (pntRes.Status == PromptStatus.OK)
            { anchorPoint = pntRes.Value; }
            else
                return;

            //trying to send an ent instead of a Polyline
            Polyline walkerPline = getPolyLine(tr, ed, msg);
            if (walkerPline == null)
                return;

            //request linked point
            ppo.Message = "\nSelect start point from walking pattern.";
            Point3d walkPoint = new Point3d();
            pntRes = ed.GetPoint(ppo);
            if (pntRes.Status == PromptStatus.OK)
            { walkPoint = pntRes.Value; }
            else
                return;


            //initiate the jig
            MutateJig jig = new MutateJig(anchorPline, walkerPline, anchorPoint,walkPoint);

            //get dynamic point to start this domino
            PromptResult res = ed.Drag(jig);
            if (res.Status == PromptStatus.OK)
            {
                //code to do stuff with dragpoint

                //get final data and dispose of the temp clone
                //might not want to make changes, in that case need to dispose of this 

                Matrix3d trans = jig.Transformation;

                walkerPline.UpgradeOpen();
                walkerPline.TransformBy(trans);
            }
            #region explaining jigger class

            //use these objects and Points to create a jig that will
            //rotate and move a polyline (or collection of ents) around another polyline
            //walking polyline should remain both parrallel at the touching points
            //as well as flipped outside of the pattern

            //track distance from start point on first polyline
                //  *instead of tracking curser from point (using some kind of nearest point formula)
                //  might work better to use the arrow keys to give the numbers to the formula
            //calculate point on corresponding polyline
            //  (moving CW or CCW, whichever is opposite curser movement)

            //use the second point to move/rotate second polyline (collection)
            //second point should share the first point and be parrallel
            //use the slope of the point on the second part compared to the slope of the point on the anchor pattern
            //  *the difference in slopes is the changing degree we need to rotate
            //  the pattern on its point, rotation is either that degree or that + 1/2 turn based on CW or CCW
            #endregion
            tr.Commit();
        }
    }

    private Polyline getPolyLine(Transaction tr, Editor ed, string msg)
    {
        //use a filter to only select polylines(not collections of entities just yet)
        //might need a collection of entities by the end, so that notches will rotate with them
        TypedValue[] filterList = new TypedValue[1];
        filterList.SetValue(new TypedValue(0, "LWPOLYLINE"), 0);
        SelectionFilter setFilter = new SelectionFilter(filterList);
        //PromptSelectionOptions pso = new PromptSelectionOptions();
        //pso.MessageForAdding = msg;

        PromptEntityOptions peo = new PromptEntityOptions(msg);
        //if user fails to select a polyline or cancels we'll return null and end the command
        PromptEntityResult select = ed.GetEntity(peo);
        if (select.Status == PromptStatus.OK)
        {
            DBObject dbo = tr.GetObject(select.ObjectId, OpenMode.ForRead) as DBObject;
            if (dbo is Polyline)
            {
                Polyline pLine = dbo as Polyline;
                return pLine;
            }
            else
            { return null; }
        }
        else
        { return null; }
    }
}

}
