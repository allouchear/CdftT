#ifndef CDFTT_GRIDPOINTS_H_INCLUDED
#define CDFTT_GRIDPOINTS_H_INCLUDED

#include<vector>

class GridPoints
{
    private:
        std::vector<int> _Lebedev_Npts;
        std::vector<int> _Lebedev_Lmax;
        std::vector<int> _Lebedev_L2max;
        int _Npts;
        int _Lmax;
        int _L2max;
        std::vector<std::vector<double>> _LebedevGridPoints;
    public:

            //! A default constructor.
            /*! This constructor is used to set all of the parameters for one LCAO on 0 or "None" value. */
        GridPoints();

            //! A real constructor.
            /*! This constructor is used to initialize the good grid. */
        GridPoints(int);

            //! A default desctructor.
            /*! We don't use it. */
        ~GridPoints() {}

            //! A normal member taking no arguments and returning a vector<int> value.
            /*! \return The Lmax of a grid. */
        std::vector<int> Lebedev_Lmax() {return _Lebedev_Lmax;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of points which compose the grid. */
        int Npts() {return _Npts;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The Lmax for the grid. */
        int Lmax() {return _Lmax;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The L2max the grid. */
        int L2max() {return _L2max;}

            //! A normal member taking no arguments and returning a vector<vector<double>> value.
            /*! \return All of the coordonates of points which compose the grid. */
        const std::vector<std::vector<double>>& LebedevGridPoints() const {return _LebedevGridPoints;}

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 6 points */
        void GridPoints6();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 14 points */
        void GridPoints14();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 26 points */
        void GridPoints26();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 38 points */
        void GridPoints38();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 50 points */
        void GridPoints50();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 74 points */
        void GridPoints74();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 86 points */
        void GridPoints86();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 110 points */
        void GridPoints110();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 146 points */
        void GridPoints146();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 170 points */
        void GridPoints170();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 194 points */
        void GridPoints194();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 230 points */
        void GridPoints230();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 266 points */
        void GridPoints266();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 302 points */
        void GridPoints302();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 350 points */
        void GridPoints350();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 434 points */
        void GridPoints434();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 590 points */
        void GridPoints590();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 770 points */
        void GridPoints770();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 974 points */
        void GridPoints974();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 1202 points */
        void GridPoints1202();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 1454 points */
        void GridPoints1454();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 1730 points */
        void GridPoints1730();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 2030 points */
        void GridPoints2030();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 2354 points */
        void GridPoints2354();

            //! A normal member taking no arguments and returning a void value.
            /*! Initialize the grid with 5810 points */
        void GridPoints5810();
};

#endif