"""ThursdayChecker
Script to help with Thursday cecking of annotated TrEMBL entries

Usage: ttc.py <file>

If <file> contains checking markup this is extracted, formatted and passed
as message body to Thunderbird.

Checking markup:
---------------

line type: ##

Example:  ## This is where the comments go
          ##  This is where the comments go
          ##   This is where the comments go

Similar to other lines types in UniProt entries, like DE, GN etc.
Multiline comments are possible but each line has to start with ##.

context specifier: #[0-99], optional, default:1, max:9

Example:  ## #5 This where the comments go
         Comment plus 5 preceding lines from entry appear in report.

          ##   #9 This where the comments go (multiline)
          ##   This where the comments go
          Comment (multiline) plus 9 preceding lines from entry appear
          in report.

After the ## line type, the context specifier allows specifying the number of
lines from the TrEMBL entries preceding the ## to be included in the final report.
When no context specifier is given, the context default to 1 line. Arbitrarily, the
maximum number of lines for context is set to 9 although double-digit figures are
parsed but ignored. If more context is specified than lines have been parsed from
current entry (e.g. 9 lines for a comment on the ID line), all lines that have been
parsed wil be returned.

Setup for Crisp:
extension.cr: extended with macro 'tcc' which saves current buffer in tcc.tmp
    location: C:\\UniProt\crisp\macros\sprot\extension.cr(cm)
tcc.bat: called by macro 'tcc', runs tcc.py on ttc.tmp
    location: C:\\UniProt\binamos\tcc.bat
"""
from collections import deque
from io import StringIO
from itertools import count
from textwrap import TextWrapper
from datetime import date
import string
import re
import sys
import os
import subprocess

#TODO: make the script handle logifles, too
class Feedback(object):
    """(Multiline) feedback plus its context from the entry."""

    _ids = count(0) #TODO: instance counter, not yet used

    def __init__(self, feedback=None, context=None, ac=None):
        """Initialize with default values."""
        self.id = next(self._ids) #TODO: not yet used
        self.feedback = StringIO(feedback) #we just write to for multiline comments
        self.context = StringIO(context)
        self.ac = ac #UniProt accession
        self.tag = None #TODO: allows re-use of a comment by referring to tag

    def add_feedback(self, aString):
        """Accumulates the feedback."""
        out = aString[2:].strip() + " "
        self.feedback.write(out)

    def add_context(self, aString):
        """Accumulates the context."""
        self.context.write(aString)

    def __str__(self):
        """String representation of instance """
        w = TextWrapper(width=75,
                        initial_indent="##   ",
                        subsequent_indent = "##   ")
        wrapped_text = w.fill(self.feedback.getvalue())
        return "{0}:\n{1}\n{2}".format(self.ac,
                                       self.context.getvalue(),
                                       wrapped_text)


class Report(object):
    """Collects all the feedback per entry"""

    def __init__(self):
        self.file = None # filename, stored for report
        self.scanner = deque(maxlen = 10) #Once we hit a ## line,
                                          #we still have 9 lines of context
        self.current_entry = None
        self.current_context = 1 #no. of lines, int 0-9
        self.all_accessions = []
        self.all_feedback = [] #instance of Feedback
        self.context_regex = re.compile("#[0-9][0-9]?") # ?: zero or 1 matches of preceding RE
        self.ac_regex = re.compile("[A-Z][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]") #UniProt AC

    def parse(self, filename):
        self.file = os.path.split(filename)[1] #just the actual file name
        with open(filename, "r", encoding='latin1') as f:
            fb_block = False
            for line in f:
                self.scanner.append(line.strip())
                ac_match = self.ac_regex.match(line) #we find this in logfiles
                if line.startswith("AC "): #new, sub or pep file
                    fb_block = False
                    acc = line[5:].split(";")
                    self.current_entry = acc[0] #only primary AC is of interest
                elif line.startswith("AC:"): #logfile
                    pass #In logfiles, this line is of no interest
                elif ac_match is not None:
                    if self.current_entry is None:
                        self.current_entry = ac_match.group(0)
                    else:
                        pass #we already have an AC; logfiles can have multiple ones
                elif line.startswith("##"):
                    if not fb_block:
                        fb_block = True
                        fb = self.__generate_feedback(ac=self.current_entry)
                        self.__reset_current_context()
                        self.__provide_context(line, fb)
                        fb.add_feedback(line)
                    else:
                        fb = self.__get_current_feedback()
                        fb.add_feedback(line)
                elif line.startswith("//"):
                    self.__reset()
                    fb_block = False
                else:
                    fb_block = False
            self.__retrofit_ACs()

    def __retrofit_ACs(self): #FIX: doesn't work for logfiles
        """Feedback for ID line, i.e. before the AC line has no AC assigned to it
        as it hasn't been read yet, so needs to be retrofitted.
        """
        for fb in self.all_feedback:
            if fb.ac is None:
                try:
                    fb_idx = self.all_feedback.index(fb)
                    fb.ac = self.all_feedback[fb_idx + 1].ac
                except IndexError:
                    fb.ac = self.all_accessions[-1] #only happens for very last entry


    def __generate_feedback(self, feedback=None, context=None, ac=None):
        fb = Feedback(feedback=feedback, context=context, ac=ac)
        self.all_feedback.append(fb)
        return fb

    def __get_current_feedback(self):
        return self.all_feedback[-1]

    def __check_for_context(self, aString):
        match = self.context_regex.search(aString)
        if match is not None and len(match.group(0)) < 3:
            self.current_context =  int(match.group(0)[1]) #2nd character of a 2 character match
        elif match is not None and len(match.group(0)) == 3:
            self.current_context = 9 #even if values > 9 are given, we limit it here
        else:
            pass

    def __extract_context(self):
        lines_of_context = len(self.scanner)
        end_idx = lines_of_context-1
        start_idx = end_idx - self.current_context
        d_list = list(self.scanner)
        if start_idx >= 0:
            context = "\n".join(d_list[start_idx:end_idx])
        else:
            context = "\n".join(d_list[0:end_idx])
        return context

    def __provide_context(self, line, feedback_instance):
        self.__check_for_context(line)
        feedback_instance.add_context(self.__extract_context())

    def __reset_current_context(self):
        self.current_context = 1

    def __reset(self):
        self.__reset_current_context()
        self.all_accessions.append(self.current_entry)
        self.current_entry = None
        self.scanner.clear()

    def _get_default_mail_client(self):
        """Query Windows registry for default client"""
        import winreg as wr
        aReg = wr.ConnectRegistry(None, wr.HKEY_CLASSES_ROOT)
        try:
            targ = r"mailto\shell\open\command"
            aKey = wr.OpenKey(aReg, targ)
            mail_client = wr.EnumValue(aKey, 0)[1]#tuple with the value at index 1
            idx = string.index(mail_client, ".exe")
            return str(mail_client[1:idx+4]) #increasing the index to inlcude .exe, no preceding/trailing "
        finally:
            wr.CloseKey(aKey)
        wr.CloseKey(aReg)

    def __send_email(self, report=None):
        """This geared towards Thunderbird as command line options
        for other client are unknown.

        """
        default_text = "Please paste report here."
        datum = date.today().strftime("%d/%m/%Y")
#        args = [self._get_default_mail_client(), "-compose"]
        args = ['C:/Program Files (x86)/Mozilla Thunderbird/thunderbird', "-compose"]
        if report is not None:
            try:
                text = report.getvalue()
            except:
                text = default_text
        else:
            text = default_text
        compose_args = "subject='Checking {0}',body='{1}'".format(datum, text)
        args.append(compose_args)
        subprocess.call(args)

    def generate_report(self, target="stdout"):
        """Summarizes all the feedback.

        Parameter:
        ----------
        target: optional
                values: stdout, email, file
                default: stdout
        """

        if target == "stdout":
            print("{0}\n\n".format(self.file))
            for fb in self.all_feedback:
                print(fb)
                print("-----")
        elif target == "email":
            report = StringIO()
            report.write("{0}:\n\n".format(self.file))
            for fb in self.all_feedback:
                report.write(fb.__str__())
                report.write("\n-----\n")
            self.__send_email(report)
        elif target == "file": #TODO: to be implemented
            pass


if __name__ == '__main__':
    if len(sys.argv) > 1:
        #target = os.getcwd() + os.sep + sys.argv[1]
        target = sys.argv[1]
    else:
        print(__doc__)
        filename = input('Which file would you like to parse? (x to quit) ')
        if filename == "x":
            sys.exit("Bye")
        else:
#            target = os.getcwd() + os.sep + filename
            target = filename
    rep = Report()
    #rep.parse("tcc_test.txt")
    rep.parse(target)
    rep.generate_report(target="email")
