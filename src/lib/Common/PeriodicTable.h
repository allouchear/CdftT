#ifndef CDFTT_PERIODICTABLE_H_INCLUDED
#define CDFTT_PERIODICTABLE_H_INCLUDED

#include <string>
#include <vector>

#include "../Common/Element.h"


    //! A periodic table class.
    /*! This class will be used in the program to have access of all the properties of an element. */
class PeriodicTable
{
    private:
        std::vector<Element> _periodic_table;
    public:

            //! A default constructor.
            /*! This create a periodic table. */
        PeriodicTable();

            //! A default desctructor.
            /*! We don't use it. */
        ~PeriodicTable();

            //! A normal member taking one argument and returning an element value.
            /*! 
                \param i an integer argument (the atomic number).
                 \return The element corresponding to the atomic number. */
        Element element(const int i) const;

            //! A normal member taking one argument and returning an element value.
            /*! 
                \param symbolOrName a std::string argument (the symbol or the name).
                 \return The element corresponding to the symbol or the name. */
        Element element(const std::string& symbolOrName) const;

            //! A normal member taking one argument.
            /*! Add an element.
                \param E an element argument. */
        void add_element(const Element& e);

            //! A normal member taking no arguments.
            /*! Add all elements and their isotopes  in the periodic table. */
        void add_all_element();

            //! A normal member taking one argument.
            /*! Add an isotope in an element. 
                \param I an isotope argument. */
        void add_isotope(const Isotope& isotope);

            //! A normal member taking no arguments.
            /*! Add all isotopes for each elements.
                \sa _add_all_element. */
        void add_all_isotope();
};

#endif
