#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
#include <fstream>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    // We assume
    //  - sequence contains at least one character
    //  - all characters in sequence are A, C, T, G, or N
    string       m_name;
    string       m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
: m_name(nm), m_sequence(sequence)
{}

// Loads a file into the appropriate genome objects. Returns true if
// load was successful, false if there was improper file formatting
bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
    // Empty the Genome vector parameter, genomes
    genomes.clear();
    
    // IMPROPER FILE FORMAT CHECKED
    //  - not starting with a name line
    //  - line starting with '>' but containing no other characters
    //  - non-name lines contain anything other than upper/lower A C T N G
    //  - no base lines after a name line
    //  - empty lines present
    
    string line;
    while (getline(genomeSource, line))
    {
        // We check to see that there is a first line
        // If so, we ensure that the first char is a '>'
        if (line.size() == 0)
            return false;
        else if (line[0] != '>')
            return false;
        
        // First char must be a '>', so we check to see if
        // the rest of the line is a name
        string name = line.substr(1);
        if (name.size() == 0)
            return false;
        
        // If the line after the name is empty, return false
        if (genomeSource.peek() == '\n')
            return false;
        
        // We have a name, so now check for a sequence
        // by looping through characters until we hit a '>'
        // or end of file
        char base;
        string sequence;
        while (genomeSource.get(base))
        {
            if (genomeSource.peek() == '>')
                break;
            
            base = toupper(base);
        
            // If there is an empty line present, we return false
            if (base == '\n')
            {
                if (genomeSource.peek() == '\n')
                    return false;
                else
                    continue;
            }
            
            // If there is a character other than upper or lower ACGTN, we return false
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T' && base != 'N')
                return false;
            
            sequence += base;
        }
        
        // If there are no base lines after a name line, we return false
        if (sequence.size() == 0)
            return false;
        
        
        // name is set to a genome's name
        // sequence is the sequence of bases
        Genome g(name, sequence);
        
        // Add the new genome to our library
        genomes.push_back(g);
    }
    
    return true;
}

int GenomeImpl::length() const
{
    return m_sequence.length();
}

string GenomeImpl::name() const
{
    return m_name;
}

// Allows user to extract a portion of the genome's sequence
// starting at a position, and running a certain length.
// Returns true if extraction was successful, false otherwise
bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    if (position < 0 || position >= m_sequence.length())
        return false;
    
    // If fragment would over-reach the available sequence
    // return false
    if (position + length > m_sequence.length())
        return false;
    
    fragment = m_sequence.substr(position, length);
    
    return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
