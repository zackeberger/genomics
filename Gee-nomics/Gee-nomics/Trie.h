//
//  GenomeMatcher.hpp
//  Gee-nomics
//
//  Created by Zack Berger on 3/6/19.
//  Copyright © 2019 Berger Inc. All rights reserved.
//

#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
using namespace std;

// Assume that characters are 1 byte and
// can be converted to integers
const int POSSIBLE_CHARACTERS = 256;

template<typename ValueType>
class Trie
{
public:
    Trie();
    // Initializes one new Node into the Trie structure and
    // points m_root at it
    
    ~Trie();
    // Frees all memory used up by the Trie
    
    void reset();
    // Frees all memory used by the Trie then
    // initializes a new empty node that m_root points to
    
    void insert(const std::string& key, const ValueType& value);
    // Associates the specified key with a specified value
    
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;
    // Searches for the values associated with a given key.
    // User can specify if values exactly matching the search key should be
    // returned (exactMatchOnly == true), or values that
    //     - match the first character of the key exactly
    //     - have a single mismatching character of the key
    //       anywhere past the first character
    
    // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    
    // Node representation:
    //  - Vector of values present in the node
    //  - Array of node pointers to children, one
    //    pointer for each possible character
    struct Node
    {
        Node()
        {
            for (int i = 0; i < POSSIBLE_CHARACTERS; i++)
                m_children[i] = nullptr;
        }
        
        vector<ValueType>   m_values;
        Node*               m_children[POSSIBLE_CHARACTERS];
    };
    
    Node*   m_root;
    
    // PRIVATE HELPER FUNCTIONS
    void freeAllNodes(Node* root);
    void findHelper(const std::string& key, bool exactMatchOnly, Node* t, vector<ValueType>& v) const;
};

template<typename ValueType>
Trie<ValueType>::Trie()
{
    m_root = new Node();
}

template<typename ValueType>
Trie<ValueType>::~Trie()
{
    freeAllNodes(m_root);
}

template<typename ValueType>
void Trie<ValueType>::reset()
{
    freeAllNodes(m_root);
    m_root = new Node();
}

template<typename ValueType>
void Trie<ValueType>::freeAllNodes(Node* root)
{
    if (root != nullptr)
    {
        for (int i = 0; i < POSSIBLE_CHARACTERS; i++)
            freeAllNodes(root->m_children[i]);

        delete root;
    }
}


template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
    string sequence = key;
    Node* curNode = m_root;         // pointer to the currrent Node
    
    while (sequence != "")
    {
        char firstChar = sequence[0];       // the first character in the current sequence
        
        // If the current node does not have a path to the current character
        // Initialize a new node that the character slot points to
        if (curNode->m_children[firstChar] == nullptr)
            curNode->m_children[firstChar] = new Node();
        
        // Make pointer to curNode move forward in the path
        curNode = curNode->m_children[firstChar];
        
        sequence = sequence.substr(1);
    }
    
    // curNode points to the node that corresponds to the key.
    // Insert the value at that node
    curNode->m_values.push_back(value);
}


template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
    vector<ValueType> temp;
    
    if (key.size() == 0)
    {
        // Return whatever is stored in the first node in the Trie, the "" node
        temp.insert(temp.end(), m_root->m_values.begin(), m_root->m_values.end());
        return temp;
    }
    
    char firstChar = key[0];
    
    // Key has not been stored in the Trie yet
    if (m_root->m_children[firstChar] == nullptr)
        return temp;
    
    // We checked that the first character matched exactly to the search key
    // so now run our findHelper function, which permits up to one mismatched character
    findHelper(key.substr(1), exactMatchOnly, m_root->m_children[firstChar], temp);
    
    return temp;
}


template<typename ValueType>
void Trie<ValueType>::findHelper(const std::string& key, bool exactMatchOnly,
                                                   Node* t, vector<ValueType>& v) const
{
    if (t == nullptr)
    {
        // Search term is not present, so do nothing and return the empty vector
        return;
    }
    
    if (key.size() == 0)
    {
        // Found a Node that a valid key maps to, so store the node's values
        // in our vector temp
        v.insert(v.end(), t->m_values.begin(), t->m_values.end());
        return;
    }
    
    char firstChar = key[0];
    
    if (exactMatchOnly)
    {
        findHelper(key.substr(1), true, t->m_children[firstChar], v);
    }
    else
    {
        // We allow up to one mismatch. Recursively call findHelper(),
        // traversing down the firstChar's path. Continue to allow up to one
        // mismatch for that call
        findHelper(key.substr(1), false, t->m_children[firstChar], v);
        
        for (int i = 0; i < POSSIBLE_CHARACTERS; i++)
        {
            // Recursively call findHelper for each character. If that path
            // continues, do not allow any more mismatches, because this was the
            // one permitted mismatch
            if (i != firstChar)
            {
                findHelper(key.substr(1), true, t->m_children[i], v);
            }
        }
    }
}


#endif // TRIE_INCLUDED
