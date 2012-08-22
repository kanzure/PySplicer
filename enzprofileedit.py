#!/usr/bin/env python3
# Takes an enzyme name and a target profile, and adds/removes specified enzyme from profile.
#  Master list is in ./lib/allenzymes.json
# ==================================================
# Program Directory implementation:
# ./codonopt.py
# ./chenzprofile.py
# ./lib/codontables/default_ecoli.json
# ./lib/codontables/default_bsubtilis.json
# ./lib/codontables/default_generic_eubacterial.json
# ./lib/allenzymes.json
# ./lib/enzymeprofiles/biobrick.json

import sys, os, json, argparse

# ================================================================================= #
ArgParser = argparse.ArgumentParser(description=('Takes comma-delimited lowercase enzyme names and a target profile, and adds/removes specified enzyme from profile.'
                                                ' Depends on "allenzymes.json", normally found in the included ./lib/ directory.'
                                                ' Profiles are created and edited in the ./lib/enzymeprofiles/ directory.'))
ArgParser.add_argument('Profile', type=str,
                       help='The name of the target profile (without .json extension). If it does not exist it will be created.')
ArgParser.add_argument('enzymes', nargs='?', type=str,
                       help='A comma-delimited list of enzyme names, preferably lowercase e.g. "ecori, bsumi".')
ArgParser.add_argument('-a', '--add', action='store_true',
                       help='Add specified enzymes to the target profile.')
ArgParser.add_argument('-r', '--remove', action='store_true',
                       help='Remove specified enzymes to the target profile.')
ArgParser.add_argument('-c', '--comment', type=str,
                       help='An optional comment to add to the profile, i.e. "For use with B.subtilis designs".')
ArgParser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose mode: creates more non-error messages.')
Args = vars(ArgParser.parse_args())
# ================================================================================= #
def verboseMsg(Msg): # Optional print statement.
    if Args['verbose']: print(Msg)

class Profile:
    def __init__(self, Name, Comment=''):
        self.scriptDir = sys.path[0] # Gives script directory. Default/imported codon tables will be here.
        self.libDir = os.path.join(self.scriptDir, 'lib')
        self.profilesDir = os.path.join(self.libDir, 'enzymeprofiles')
        self.workingDir = os.getcwd() # Gives current working directory.
        self.name = Name
        self.profilepath = os.path.join(self.profilesDir, self.name+'.json')
        try:
            self.profile = self.getProfile() # Either loads or creates if it cannot.
        except AssertionError:
            print("Error: Designated profile not correctly formatted: must have string-type 'name' entry and list-type 'EnzymeList' entry.")
            sys.exit(1)
        try:
            with open(os.path.join(self.libDir, 'allenzymes.json'), mode='r') as MasterFile:
                self.masterdict = json.loads(MasterFile.read())
        except IOError:
            print("Unable to find 'allenzymes.json' in lib directory. Aborting.")
            sys.exit(1)
        if Comment: # Empty strings eval false, so by default this does nothing.
            self.comment = Comment
            self.profile['comment'] = self.comment
        elif 'comment' not in self.profile.keys():
            self.profile['comment'] = '' #Just to avoid keyerror bugs if we ever call comment for some reason..

    def addEnzyme(self,Enzyme):
        def AddFromMaster(AEnzyme):
            if AEnzyme in self.masterdict['EnzymeList']:
                self.profile[AEnzyme] = self.masterdict[AEnzyme]
                return True
            else:
                print("Error: Unable to add enzyme '"+Enzyme+"' to profile: Enzyme not found. Continuing.")
                return False
        if Enzyme in self.profile['EnzymeList']:
            if Enzyme not in self.profile.keys():
                verboseMsg("Enzyme found in profile's 'EnzymeList' but is absent; adding enzyme entry.")
                AddFromMaster(Enzyme)
            verboseMsg('Enzyme already in profile; skipping.')
            return None # Abort
        else:
            if AddFromMaster(Enzyme):
                self.profile['EnzymeList'].append(Enzyme)
                verboseMsg("Added enzyme '" + Enzyme + "' to profile.")

    def removeEnzyme(self,Enzyme):
        def RemoveListEntryByName(Entry, List):
            EntryIndex = List.index(Entry)
            del(List[EntryIndex])
        if Enzyme in self.profile['EnzymeList']:
            try:
                del(self.profile[Enzyme])
            except KeyError:
                verboseMsg('Enzyme found in "EnzymeList" but not in Profile. Removing entry.')
            RemoveListEntryByName(Enzyme, self.profile['EnzymeList'])
            verboseMsg("Enzyme '"+Enzyme+"' deleted from Profile.")
            if Enzyme in self.profile['EnzymeList']:
                verboseMsg("Double entry detected for "+ Enzyme +", removing second entry. Others may remain.")
                RemoveListEntryByName(Enzyme, self.profile['EnzymeList'])
        else:
            verboseMsg("Enzyme '"+ Enzyme +"' not found, not deleted. Continuing.")

    def addEnzymeList(self,EnzymeList):
        assert isinstance(EnzymeList, list)
        for Enzyme in EnzymeList:
            assert isinstance(Enzyme, str)
            self.addEnzyme(Enzyme.lower())

    def removeEnzymeList(self,EnzymeList):
        assert isinstance(EnzymeList, list)
        for Enzyme in EnzymeList:
            assert isinstance(Enzyme, str)
            self.removeEnzyme(Enzyme.lower())

    def createProfile(self):
        return {'name'		: self.name,
                'EnzymeList'	: []}

    def getProfile(self):
        try:
            with open(self.profilepath, mode='r') as ProfileFile:
                Profile = json.loads(ProfileFile.read())
        except IOError:
            verboseMsg("Couldn't find specified profile in enzymeprofiles directory: creating new profile.")
            Profile = self.createProfile()
        assert 'EnzymeList' in Profile.keys()
        assert isinstance(Profile['EnzymeList'], list)
        assert 'name' in Profile.keys()
        assert isinstance(Profile['name'], str)
        return Profile

    def saveProfile(self):
        try:
            with open(self.profilepath, encoding='utf-8', mode='w') as ProfileFile:
                ProfileFile.write(json.dumps(self.profile, indent=1))
        except:
            print("Unable to save to ", self.profilepath, " - Does path exist? Do you have write permissions?")
            sys.exit(1)

# Expect args to contain: Profile (str), enzymes (comma-delimited str), add (True/None), remove (True/None), verbose (True/None)
if Args['add'] and Args['remove']:
    print("Error: You cannot specify the '-a' and '-r' flags simultaneously!")
    sys.exit(1)
if Args['comment']:
    CommentArg = Args['comment']
else:
    CommentArg = 'No Comment'
NewProfile = Profile(Args['Profile'], CommentArg)

Enzymes = []
for Enzyme in Args['enzymes'].split(","):
    Enzymes.append(Enzyme.strip().lower())

if Args['add']:
    NewProfile.addEnzymeList(Enzymes)
    NewProfile.saveProfile()
elif Args['remove']:
    NewProfile.removeEnzymeList(Enzymes)
    NewProfile.saveProfile()
else:
    print("You must specify either '-a' or '-r' to create or modify profiles.")
    sys.exit()

