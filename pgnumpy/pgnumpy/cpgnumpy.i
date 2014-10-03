%module cpgnumpy
%{
#include "cpgnumpy.h"
%}

%typemap(throws) const char * %{
    PyErr_SetString(PyExc_RuntimeError, $1);
    SWIG_fail;
%}

class cPgNumpy {

    public:

        cPgNumpy() throw (const char*);
        cPgNumpy(string conninfo) throw (const char*);
        cPgNumpy(const char* conninfo) throw (const char*);
        ~cPgNumpy();

        void open() throw (const char*);
        void open(string conninfo) throw (const char*);
        void open(const char* conninfo) throw (const char*);


        // Execute a query
        void execute(string query_string) throw (const char*);
        void execute(const char* query_string) throw (const char*);

        // fetch results as an array
        PyObject* fetchall() throw (const char*);


        // Write the results to a file or stdout. Requires running
        // use_text() to get results in text format
        long long write(const char* filename=NULL) throw (const char*);

        // return results as text instead of binary.  Note, currently we don't
        // support fetching the results into numpy arrays for text retrevial,
        // but setting text is required for the write() method to write to a
        // file or stdout
        void use_text();

        // should we use array decorators for array writing? This will
        // use { } to delineate array dimensions.  This is required
        // for COPYing back into a table.  DEFAULT is true.
        void set_decorators(bool decorate);



        // Execute the query and use a cursor to retrieve the results in
        // chunks, writing to a file or by default stdout.  Unlike write(),
        // does *not* require use_text() to be set.
        //
        // This does not use the query execution infrastructure of cPgNumpy,
        // but rather is rolled specially for the cursor
        //
        // The return value is the number of rows written
        long long execwrite(
                const char* query,
                const char* filename=NULL) throw (const char*);
        // size of fetch chunks when doing an execwrite
        // DEFAULT is 1000
        void set_fetch_count(long long fetchcount);


        long long ntuples() {
            return mNtuples;
        }
        long long nfields() {
            return mNfields;
        }
        int status() {
            return mResultStatus;
        }

        // explicitly set fields lengths.  This is useful for varchar types,
        // since for those types all the returned values must be scanned to get
        // the *largest*, and the output array field will have those
        // dimensions.  With this setting you can shortcut that.
        void set_field_lengths(PyObject* flenDict) throw (const char*);
        void clear_field_lengths() throw (const char*);

        void clear();
        void close();



};
