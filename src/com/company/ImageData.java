package com.company;
import org.ejml.data.DMatrixRMaj;

import java.awt.*;
import java.util.ArrayList;

import static org.ejml.dense.row.CommonOps_DDRM.mult;

public class ImageData
{
    public static class ImageData2D
    {
        public double width;
        public double height;
        public ArrayList<PairPoint2D> pair_points;
    }

    public static class Point2D
    {
        public double x;
        public double y;
    }

    public static class PairPoint2D
    {
        public Point2D image_point;
        public Point2D model_point;
    }

    public static class Point3D
    {
        public double x;
        public double y;
        public double z;
    }

    public static class PairPoint3D
    {
        public Point3D image_point;
        public Point3D model_point;
    }

    public static class ImageData3D
    {
        public double width;
        public double height;
        public ArrayList<PairPoint3D> pair_points;
    }

    public static Point3D point_2d_to_point_3d (Point2D p_in)
    {
        Point3D p = new Point3D();
        p.x = p_in.x;
        p.y = p_in.y;
        p.z = 0;
        return p;
    }

    public static class Point4D
    {
        public double x;
        public double y;
        public double z;
        public double w;
    }

    public static DMatrixRMaj point_2d_to_colum_vector (Point2D p)
    {
        DMatrixRMaj mat = new DMatrixRMaj(2, 1);
        mat.set(0, 0, p.x);
        mat.set(1, 0, p.y);
        return mat;
    }

    public static Point2D colum_vector_2d_to_point_2d (DMatrixRMaj m)
    {
        if (m.getNumCols() != 1 && m.getNumRows() != 2)
        {
            throw new AssertionError("Matrix dimension does not fit");
        }
        Point2D out = new Point2D();
        out.x = m.get(0,0);
        out.y = m.get(1,0);
        return out;
    }

    public static DMatrixRMaj point_3d_to_colum_vector (Point3D p)
    {
        DMatrixRMaj mat = new DMatrixRMaj(3, 1);
        mat.set(0, 0, p.x);
        mat.set(1, 0, p.y);
        mat.set(2, 0, p.z);
        return mat;
    }

    public static Point3D colum_vector_3d_to_point_3d (DMatrixRMaj m)
    {
        if (m.getNumCols() != 1 && m.getNumRows() != 3)
        {
            throw new AssertionError("Matrix dimension does not fit");
        }
        Point3D out = new Point3D();
        out.x = m.get(0,0);
        out.y = m.get(1,0);
        out.z = m.get(2,0);
        return out;
    }

    public static DMatrixRMaj point_4d_to_colum_vector (Point4D p)
    {
        DMatrixRMaj mat = new DMatrixRMaj(3, 1);
        mat.set(0, 0, p.x);
        mat.set(1, 0, p.y);
        mat.set(2, 0, p.z);
        mat.set(3, 0, p.w);
        return mat;
    }

    public static Point4D colum_vector_4d_to_point_4d (DMatrixRMaj m)
    {
        if (m.getNumCols() != 1 && m.getNumRows() != 4)
        {
            throw new AssertionError("Matrix dimension does not fit");
        }
        Point4D out = new Point4D();
        out.x = m.get(0,0);
        out.y = m.get(1,0);
        out.z = m.get(2,0);
        out.w = m.get(3,0);
        return out;
    }

    public static DMatrixRMaj build_normalization_matrix(double w, double h)
    {
        DMatrixRMaj n = new DMatrixRMaj(3,3);
        n.set(0, 0, 2.0/ w);
        n.set(1, 1, 2.0 / h);
        n.set(0, 2, -1.0);
        n.set(1,2 , -1.0);
        return n;
    }

    public static DMatrixRMaj normalize_point(Point2D point, DMatrixRMaj normalization_matrix )
    {
        Point3D point_3d = point_2d_to_point_3d(point);
        DMatrixRMaj point_3d_as_mat = point_3d_to_colum_vector(point_3d);
        DMatrixRMaj out = new DMatrixRMaj(normalization_matrix.getNumRows(), 3);
        mult(normalization_matrix, point_3d_as_mat, out);
        return out;
    }

    public static ImageData3D normalize_data_set(ImageData2D data_set)  {
        DMatrixRMaj normalization_matrix = build_normalization_matrix(data_set.width, data_set.height);
        ArrayList<PairPoint3D> normalized_points  = new ArrayList<>();

        for (int i = 0; i < data_set.pair_points.size(); i++)
        {
            PairPoint2D p = data_set.pair_points.get(i);
            Point3D new_model_point = point_2d_to_point_3d(p.model_point);
            Point3D new_image_point = colum_vector_3d_to_point_3d(normalize_point(p.image_point, normalization_matrix));
            PairPoint3D new_pair = new PairPoint3D();
            new_pair.model_point = new_model_point;
            new_pair.image_point = new_image_point;

            normalized_points.add(new_pair);
        }

        ImageData3D r = new ImageData3D();
        r.width = data_set.width;
        r.height = data_set.height;
        r.pair_points = normalized_points;
        return r;
    }

}
