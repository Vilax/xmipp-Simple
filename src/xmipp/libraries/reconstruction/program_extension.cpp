/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosam@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "program_extension.h"

//Needed includes for instantiate programs
#include <data/filters.h>


#include <data/mask.h>
#include <classification/analyze_cluster.h>

void runSystem(const String &program, const String &arguments, bool useSystem) {
	if (useSystem) {
		String cmd = formatString("%s %s", program.c_str(), arguments.c_str());
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
	} else {
		runProgram(program, arguments);
	}
}

int runProgram(XmippProgram * program, const String &arguments, bool destroy)
{
    if (program == NULL)
        REPORT_ERROR(ERR_PARAM_INCORRECT, "Received a NULL as program pointer");
    program->read(arguments);
    int retCode = program->tryRun();
    if (destroy)
        delete program;
    return retCode;
}

int runProgram(const String &programName, const String &arguments)
{
    XmippProgram * program = getProgramByName(programName);
    return runProgram(program, arguments);
}

XmippProgram * getProgramByName(const String &programName)
{
    if (programName == "xmipp_mask")
        return new ProgMask();


    if (programName == "xmipp_classify_analyze_cluster")
        return new ProgAnalyzeCluster();


    return NULL;
}


