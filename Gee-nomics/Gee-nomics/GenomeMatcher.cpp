//
//  GenomeMatcher.cpp
//  Gee-nomics
//
//  Created by Zack Berger on 3/6/19.
//  Copyright © 2019 Berger Inc. All rights reserved.
//

#include "provided.h"
#include "Trie.h"
#include <string>
#include <utility>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
using namespace std;

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    int                           m_minSearchLength;
    vector<Genome>                m_genomeLibrary;
    Trie<pair<int, int>>          m_sequencedDNA;
};

bool genomeMatchCompare(const GenomeMatch& a, const GenomeMatch& b);

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
: m_minSearchLength(minSearchLength)
{}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    m_genomeLibrary.push_back(genome);

    // We sequence the genome's DNA by the minimum search length
    // and stop when there is no more of the minSearchLength bases
    // left to extract
    for (int i = 0; ; i++)
    {
        string substring = "error";
        genome.extract(i, m_minSearchLength, substring);
        
        if (substring == "error")
            break;
        
        // substring was successfully set, so we create a value for that sequence
        // representing the genome in which, and position at, the sequence was found.
        // We start counting at genome 1, 2, 3 and so forth
        pair<int, int> genomeANDposition(m_genomeLibrary.size(), i);
        
        // Insert the search key seqeunce and pair into our Trie
        m_sequencedDNA.insert(substring, genomeANDposition);
    }
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minSearchLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    if (fragment.size() < minimumLength)
        return false;
    
    if (minimumLength < m_minSearchLength)
        return false;
    
    matches.clear();
    
    string fragSearchKey = fragment.substr(0, m_minSearchLength);
    
    // matchLocations has the genome # and positions of each spot
    // that the minSearchLength (as stored by the GenomeMatcher
    // object) amount of the fragment was found
    
    // each pair holds the genome number (indexed by the library), and position
    vector<pair<int,int>> matchLocations = m_sequencedDNA.find(fragSearchKey, exactMatchOnly);
    
    // No matches between fragment and any
    // segment of any genome in the library
    if (matchLocations.size() == 0)
        return false;
    
    // Declare unordered_map that will store the best DNA
    // match for each genome
    unordered_map<int, DNAMatch> match;
    
    // For each match location...
    for (int i = 0; i < matchLocations.size(); i++)
    {
        // snipped indicates if we should allow a character mismatch when building up
        // the piece of DNA
        bool snipped = exactMatchOnly;
        
        int currentGenome = matchLocations[i].first;
        int currentPosition = matchLocations[i].second;
        
        int actualLength = m_minSearchLength;
        
        // Start with the piece of DNA with length minSearchLength, stored in the
        // matchLocations vector. We will now build up the piece, base by base,
        // checking it against the fragment and seeing if it is valid.
        for ( ; actualLength <= fragment.size(); actualLength++)
        {
            string substring = "error";
            m_genomeLibrary[currentGenome - 1].extract(currentPosition, actualLength, substring);
            // The currentGenome - 1 correction is because we started counting at
            // genome 1, 2, 3, and so forth when we added genomes.
            
            if (substring == "error")
                break;
            
            // substring is set to the DNA piece, with an extra base
            if (substring == fragment.substr(0, actualLength))
            {
                if (actualLength == fragment.size())
                    break;
                else
                    continue;
            }
            else
            {
                // If we permit one character mismatch
                if (!snipped)
                {
                    if (substring.size() == fragment.size())
                        break;
                    
                    snipped = true;
                    continue;
                }
                else        // No more character mismatches allowed
                {
                    if (substring[actualLength - 1] != fragment.at(actualLength - 1))   // Character mismatch!!!
                    {
                        actualLength--;
                        break;
                    }
                    
                    if (substring.size() == fragment.size())
                        break;
                }
            }
        }

        // actualLength is the length of the DNA piece that matches the fragment
        if (actualLength < minimumLength)
            continue;
        
        DNAMatch d;
        d.genomeName = m_genomeLibrary[currentGenome - 1].name();
        d.length = actualLength;
        d.position = currentPosition;
        
    
        // For each genome, we will store the piece with the best match in the match map.
        unordered_map<int, DNAMatch>::iterator it = match.find(currentGenome);
        
        if (it != match.end())
        {
            if (match[currentGenome].length < d.length)
            {
                match[currentGenome] = d;
            }
        }
        else
            match.insert(make_pair(currentGenome, d));
    }
    
    // Transfer our collection of DNA matches to the matches vector
    unordered_map<int, DNAMatch>::iterator it = match.begin();
    for (; it != match.end(); it++)
    {
        matches.push_back( (it->second) );
    }
    
    if (matches.size() == 0)
        return false;
    
    return true;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    if (fragmentMatchLength < m_minSearchLength)
        return false;
    
    results.clear();
    unordered_map<string, int> genomeMatches;       // Maps a genome name to the number of times a sequence occurred in it
    
    // numSequences is the number of sequences we will consider for analysis
    int numSequences = query.length() / fragmentMatchLength;
    
    if (numSequences == 0)
        return false;

    // Loop through each sequence of DNA
    for (int i = 0; i < numSequences; i++)
    {
        string sequence = "error";
        query.extract( (i * fragmentMatchLength), fragmentMatchLength, sequence);

        if (sequence == "error")
            break;
        
        vector<DNAMatch> matches;
        findGenomesWithThisDNA(sequence, fragmentMatchLength, exactMatchOnly, matches);
        
        // For each sequence of DNA, map the number of times it was found
        // in each Genome
        for (int i = 0; i < matches.size(); i++)
        {
            genomeMatches[matches[i].genomeName]++;
        }
    }
    
    // For each genome in the library, compute the number of matching sequences found
    // divided by the total number of sequences.
    // If its match percent exceeds the threshold, add it to the results vector
    for (int i = 0; i < m_genomeLibrary.size(); i++)
    {
        Genome cur = m_genomeLibrary[i];
        
        double p = ( genomeMatches[cur.name()] / static_cast<double>(numSequences) ) * 100;

        if (p > matchPercentThreshold)
        {
            GenomeMatch g;
            g.genomeName = cur.name();
            g.percentMatch = p;
            results.push_back(g);
        }
    }
    
    if (results.size() == 0)
        return false;
    
    // Vector should be returned in order of highest match, with ties broken
    // by the alphabetical ordering of the genome name
    sort(results.begin(), results.end(), &genomeMatchCompare);

    return true;
}

bool genomeMatchCompare(const GenomeMatch& a, const GenomeMatch& b)
{
    if (a.percentMatch > b.percentMatch)
        return true;
    else if (a.percentMatch < b.percentMatch)
        return false;
    else
    {
        if (a.genomeName < b.genomeName)
            return true;
        else
            return false;
    }
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}

