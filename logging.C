/*
	Author: John Nemeth
	Description: implementation file for logging.h
	Sources: Class material, online documentation for exception library
*/
#include <logging.h>
#include <string.h>
#include <iostream>

using std::cerr;
using std::endl;

// exception handling //////////////////////////////////////////////
DataFlowException::DataFlowException(const char *type, const char *error) {
	strcpy(msg, type);
	strcat(msg, error);
	Logger::LogEvent(msg);
}

// logging /////////////////////////////////////////////////////////

FILE * Logger::logger = fopen("logger", "w");

void Logger::LogEvent(const char * event) {
	fprintf(logger,"%s\n", event);
}

void Logger::LogEvent(const char * event, const char * from) {
	fprintf(logger, "%s: %s\n", from, event);
}

void Logger::LogEvent(const char * event, int num) {
	fprintf(logger, "%s: %d\n", event, num);
}

void Logger::Finalize() {
	fclose(Logger::logger);
}


