/*
	Original Author(s): CIS330 Teaching Staff
	Revision AUthor: John Nemeth
	Description: Header file for logging.C
	Sources: Class material and online resources for exception library
*/

#ifndef LOGGING_H
#define LOGGING_H

#include <exception>
#include <stdio.h>

using std::exception;


class DataFlowException : public exception
{
  public:
                         DataFlowException(const char *type, const char *error);
    virtual const char  *what() const throw() { return msg; };

  protected:
    char        msg[1024];
};


class Logger
{
  public:
    static void     LogEvent(const char *event);
    static void     LogEvent(const char *event, const char *from);
    static void     LogEvent(const char *event, int num);
    static void     Finalize();

  private:
    static   FILE  *logger;
};

#endif
