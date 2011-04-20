/*
 DelFEM (Finite Element Analysis)
 Copyright (C) 2009  Nobuyuki Umetani    n.umetani@gmail.com
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdio.h>
#include <iostream>
#include <stack>

#include "delfem/cad/cad_svg.h"
#include "delfem/objset.h"


void XMLSkipComment(FILE*& fp,char* buff,const unsigned int BUFF_SIZE,char*& ptr_in,char*& ptr_out)
{
  std::cout << "XMLSkipComment" << std::endl;
  assert( *ptr_in == '<' && *(ptr_in+1) == '!' && *(ptr_in+2) == '-' && *(ptr_in+3) == '-' ); 
  {
    char* tp = ptr_in;
    for(;;){
      char* tp0 = strchr(tp,'>');
//      std::cout << "hoge0 " << tp0 << " "<< tp0-1 << " " << tp0-2 << std::endl;
      if( tp0 != NULL ){
        if( *(tp0-1) == '-' && *(tp0-2) == '-' ){ ptr_out = tp0; return; }
        tp = tp0+1;
        continue;
      }
      break;
    }
  }
  while( fgets(buff,BUFF_SIZE,fp) != NULL ){
    char* tp = buff;    
    for(;;){
      char* tp0 = strchr(tp,'>');
//      std::cout << "hoge1 " << tp0 << std::endl;      
      if( tp0 != NULL ){
        if( *(tp0-1) == '-' && *(tp0-2) == '-' ){ ptr_out = tp0; return; }
        tp = tp0+1;
        continue;
      }
      break;
    }    
  }      
}



class CXMLClass
{  
public:
public:
  std::string name;
  std::vector<unsigned int> aIdMemberClass_;
  std::map<std::string,std::string> mapKeyValue_;
};

class CFileReaderXML
{
public:
  CFileReaderXML(const std::string& fname){
    is_fail = true;
    fp = fopen(fname.c_str(),"r");
    if( fp == NULL ){ return; }
    is_fail = false;
    BUFF_SIZE = 512*128;
    buff = new char [BUFF_SIZE];
    fgets(buff,BUFF_SIZE,fp);
    cp = buff;
  }
  ~CFileReaderXML(){
    if( fp != NULL ){ fclose(fp); }
    delete[] buff;
  }
  void AddKeyValue(char* chs, const char* che, std::vector< std::pair<std::string,std::string> >& aStr){
    unsigned int nch = che-chs;
    unsigned int ich_ks=0;
    for(;;){
      for(;ich_ks<nch;ich_ks++){
        if( chs[ich_ks] != ' ' ) break;
      }
      if( ich_ks == nch ) return;
      unsigned int ich_ke=ich_ks+1;
      for(;ich_ke<nch;ich_ke++){
        if( chs[ich_ke] == '=' ) break;
      }    
      if( ich_ke == nch ) return;
      chs[ich_ke] = '\0';
      std::string key(chs+ich_ks);
      chs[ich_ke] = '=';      
      ////
      unsigned int ich_vs = ich_ke;
      for(;ich_vs<nch;ich_vs++){
        if( chs[ich_vs] == '"' ) break;
      }
      if( ich_vs == nch ) return;
      ich_vs += 1;
      unsigned int ich_ve = ich_vs;
      for(;ich_ve<nch;ich_ve++){
        if( chs[ich_ve] == '"' ) break;
      }      
      chs[ich_ve] = '\0';
      std::string val(chs+ich_vs);
      chs[ich_ve] = '"';
      /////
      aStr.push_back( std::make_pair(key,val) );
      /////
      ich_ks = ich_ve+1;
    }
  }
  bool SearchChr(char chr,std::vector< std::pair<std::string,std::string> >& aStr){
    aStr.clear();
    char* tp0 = strchr(cp,chr);
    if( tp0 != NULL ){ 
      AddKeyValue(cp,tp0,aStr);
      cp = tp0; 
      return true; 
    }
    AddKeyValue(cp,cp+strlen(cp),aStr);
    while( fgets(buff,BUFF_SIZE,fp) != NULL ){
      cp = buff;
      char* tp0 = strchr(cp,chr);
      if( tp0 != NULL ){ 
        AddKeyValue(cp,tp0,aStr);      
        cp = tp0; 
        return true; 
      }
      AddKeyValue(cp,cp+strlen(cp),aStr);
    }
    return false;
  }
  bool SearchChr(char chr){
    char* tp0 = strchr(cp,chr);
    if( tp0 != NULL ){ cp = tp0; return true; }
    while( fgets(buff,BUFF_SIZE,fp) != NULL ){
      cp = buff;
      char* tp0 = strchr(cp,chr);
      if( tp0 != NULL ){ cp = tp0; return true; }      
    }
    return false;
  }  
  bool CmpStr(const std::string& str,int ioff=0){
    for(int i=ioff;i<ioff+str.size();i++){
      if( cp[i] != str[i] ) return false;
    }
    return true;
  }
  char GetChr(int ioff=0) const { 
    if( cp+ioff-buff < 0 ){ return '\0'; }  // ouf of memory
    return *(cp+ioff); 
  }
  void Move(int ioff){
    cp += ioff;
  }
  std::string GetStrThisLine(const std::string& seg){
    const unsigned int ncp = strlen(cp);
    const unsigned int nseg = seg.size();
    unsigned int icps = 0;
    for(;icps<ncp;icps++){
      bool iflgs = true;
      for(unsigned int iseg=0;iseg<nseg;iseg++){
        if( cp[icps] == seg[iseg] ) iflgs = false;
      }
      if( iflgs ) break;
    }
    unsigned int icpe = icps+1;
    for(;icpe<ncp;icpe++){
      bool iflge = false;
      for(unsigned int iseg=0;iseg<nseg;iseg++){
        if( cp[icpe] == seg[iseg] ) iflge = true;
      }
      if( iflge ) break;
    }    
    char ach[256*32];
    strncpy(ach,cp+icps,icpe-icps);
    *(ach+icpe-icps) = '\0';    
    cp = cp+icpe;
    return std::string(ach);
  }
public:
  bool is_fail;
  FILE* fp;
  unsigned int BUFF_SIZE;
  char* buff;
  char* cp;
};

/*
unsigned int ParseXML(Com::CObjSet<CXMLClass>& xml,FILE*& fp,char* buff,const unsigned int BUFF_SIZE,const char* ptr_in,char*& ptr_out){
  assert( *ptr_in == '<' ); 
  unsigned int id_this = xml.AddObj(std::make_pair(0,CXMLClass()));
  const bool is_header = ( *(ptr_in+1) == '?' );
    
  std::string name;
  { // parse this tag 
    char* tp0 = strchr(ptr_in,'>');
    if( tp0 != NULL ){ 
      ptr_out = tp0; 
      *tp0 = '\0';
      name = std::string(ptr_in+1);
      xml.GetObj(id_this).name = name;
      *tp0 = '>';
      if( is_header ){
        if( *(tp0-1) == '?' ){ return id_this; }
      }      
      if( *(tp0-1) == '/' ){
        std::cout << " new class closed tag " << name << std::endl;
        return id_this;
      }
      std::cout << " new class open tag: " << name << std::endl; 
    }
    else{
      name = std::string(ptr_in+1);
      xml.GetObj(id_this).name = name;      
      while( fgets(buff,BUFF_SIZE,fp) != NULL ){
        char* tp = buff;
        char* tp1 = strchr(tp,'>');
        if( tp1 != NULL ){           
          ptr_out = tp1; 
          if( is_header ){
            if( *(tp1-1) == '?' ){ return id_this; }
          }          
          if( *(tp1-1) == '/' ){ 
            std::cout << " new class closed tag " << name << std::endl;            
            return id_this;
          }
          else{
            std::cout << " new class open tag : " << name << std::endl;            
            break;
          }
        }          
      }
    }
  }
  
  //  ptr_out should be end of tag
  {    
    char* tp = ptr_out+1;
    for(;;){
      char* tp1 = strchr(tp,'<');
      if( tp1 == NULL ){ break; }
      if( *(tp1+1) == '/' ){
        char* tp2 = strchr(tp1,'>');
        ptr_out = tp2;
        return id_this;
      }
      if( *(tp1+1) == '!' && *(tp1+2) == '-' && *(tp1+3) == '-' ){
        char* tp_out;
        XMLSkipComment(fp,buff,BUFF_SIZE,tp1,tp_out);
        assert( *tp_out == '>' );
        tp = tp_out+1;
        continue;
      }
      char* tp2; 
      unsigned int id = ParseXML(xml,fp,buff,BUFF_SIZE,tp1,tp2);      
      std::cout << "member class in " << xml.GetObj(id_this).name << " : " << xml.GetObj(id).name << std::endl;
      assert( *tp2 == '>' );
      tp = tp2+1;
    }  
  }

  while( fgets(buff,BUFF_SIZE,fp) != NULL ){
    char* tp = buff;
    for(;;){
      char* tp1 = strchr(tp,'<');
      if( tp1 == NULL ){ break; }
      if( *(tp1+1) == '/' ){
        char* tp2 = strchr(tp1,'>');
        ptr_out = tp2;
        return id_this;
      }
      if( *(tp1+1) == '!' && *(tp1+2) == '-' && *(tp1+3) == '-' ){
        char* tp_out;
        XMLSkipComment(fp,buff,BUFF_SIZE,tp1,tp_out);
        assert( *tp_out == '>' );
        tp = tp_out+1;
        continue;
      }
      char* tp2; 
      unsigned int id = ParseXML(xml,fp,buff,BUFF_SIZE,tp1,tp2);      
      std::cout << "member class in " << xml.GetObj(id_this).name << " : " << xml.GetObj(id).name << std::endl;
      assert( *tp2 == '>' );
      tp = tp2+1;
    }
    char* tp1 = strchr(tp,'>');
    if( tp1 != NULL ){             
      ptr_out = tp1;
      assert( *(tp1-1) != '/' );
      return id_this;
      break;
    }
  }      
  return id_this;
}
*/
/*
unsigned int ReadXML(const std::string& fname, Com::CObjSet<CXMLClass>& xml)
{  
  FILE* fp;
  fp = fopen(fname.c_str(),"r");
  if( fp == NULL ){ return 0; }
  std::cout << "file was open " << fname << std::endl;  
  ////
  xml.Clear();
  const unsigned int BUFF_SIZE = 512;
  char buff[BUFF_SIZE];  
  while( fgets(buff,BUFF_SIZE,fp) != NULL ){
    char* tp_in = strchr(buff,'<');
    if( tp_in != NULL && *(tp_in+1) == '!' && *(tp_in+2) == '-' && *(tp_in+3) == '-' ){
      char* tp_out;
      XMLSkipComment(fp,buff,BUFF_SIZE,tp_in,tp_out);
      assert( *tp_out == '>' );
      tp_in = strchr(tp_out+1,'<');
    }
    if( tp_in != NULL ){
      char* tp_out;
      ParseXML(xml,fp,buff,BUFF_SIZE,tp_in,tp_out);
    }
  }  
  fclose(fp);    
}*/


unsigned int ParseXML(Com::CObjSet<CXMLClass>& xml, CFileReaderXML& fin){
  assert( fin.GetChr() == '<' ); 
  unsigned int id_this = xml.AddObj(std::make_pair(0,CXMLClass()));
  const bool is_header = ( fin.GetChr(1) == '?' );
  std::string name;
  { // extract tag name
    if( is_header ){ fin.Move(2); }
    else{ fin.Move(1); }
    name = fin.GetStrThisLine(" >\n");
    xml.GetObj(id_this).name = name;
  }
  //    std::cout << "parse xml class start : " << name << std::endl;      
  { // parse this tag 
    std::vector< std::pair<std::string,std::string> > aStr;
    if( !fin.SearchChr('>',aStr) ){ std::cout << "error : tag must be closed" << std::endl; }
    for(unsigned int istr=0;istr<aStr.size();istr++){
      std::map<std::string,std::string>& map = xml.GetObj(id_this).mapKeyValue_;
      map.insert( aStr[istr] );
    }
    if( is_header ){
      if( fin.GetChr(-1) == '?' ){ // header end here
        return id_this; 
      }
    }      
    if( fin.GetChr(-1) == '/' ){  // self-enclosing tag
      return id_this;
    }
  }
  
  { // this tag is enclosing tag
    for(;;){
      if( !fin.SearchChr('<') ){ return id_this; } // eof
      if( fin.GetChr(1) == '/' ){ // this tag ends here
        fin.Move(2);
        std::string name_e = fin.GetStrThisLine(" >\n");
        assert( name == name_e );
        fin.SearchChr('>');
        return id_this;
      }
      if( fin.CmpStr("<!--") ){ // this is comment
        fin.SearchChr('>');
        continue;
      }
      assert( fin.GetChr() == '<' );
      unsigned int id = ParseXML(xml,fin);
      xml.GetObj(id_this).aIdMemberClass_.push_back(id);
      assert( fin.GetChr() == '>' );
    }  
  }
  return id_this;
}


void ReadXML(const std::string& fname, Com::CObjSet<CXMLClass>& xml)
{  
  xml.Clear();
  CFileReaderXML fin(fname);
  for(;;){
    std::vector< std::pair<std::string,std::string> > aStr0;    
    if( !fin.SearchChr('<',aStr0) ) break;
    if( fin.CmpStr("<!--") ){
      std::vector< std::pair<std::string,std::string> > aStr1;      
      fin.SearchChr('>',aStr1);
      continue;
    }
    assert( fin.GetChr() == '<' );
    ParseXML(xml,fin);
  }
  const std::vector<unsigned int>& aId = xml.GetAry_ObjID();
  for(unsigned int iid=0;iid<aId.size();iid++){
    unsigned int id = aId[iid];
    const CXMLClass& cxml = xml.GetObj(id);
    std::cout << "id : " << id << std::endl;
    std::cout << "name : " << cxml.name << std::endl;
    std::cout << "member class : ";
    for(unsigned int isub=0;isub<cxml.aIdMemberClass_.size();isub++){
      std::cout << cxml.aIdMemberClass_[isub] << " ";
    }
    std::cout << std::endl;
    ////
    std::map<std::string,std::string>::const_iterator itr = cxml.mapKeyValue_.begin();
    std::cout << "property" << std::endl;
    for(;itr!=cxml.mapKeyValue_.end();itr++){
      std::cout << "  " << itr->first << " " << itr->second << std::endl;
    }
    std::cout << std::endl;
  }
}


std::vector<unsigned int> Cad::ReadSVG_AddLoopCad(const std::string& fname, CCadObj2D& cad_2d, double scale)
{
  std::vector<unsigned int> res;
  Com::CObjSet<CXMLClass> xml;
  ReadXML(fname,xml);
  ////
  const std::vector<unsigned int>& aId = xml.GetAry_ObjID();  
  for(unsigned int iid=0;iid<aId.size();iid++){
    unsigned int id_xml = aId[iid];
    const CXMLClass& cxml = xml.GetObj(id_xml);
    if( cxml.name != "path" ) continue;
    std::map<std::string,std::string>::const_iterator itr = cxml.mapKeyValue_.find("d");
    std::cout << itr->second << std::endl;
    std::vector<std::string> aStr;
    {
      char chpath[256*128];    
      strcpy(chpath,itr->second.c_str());      
      char* pch = strtok(chpath," ,");
      while(pch!=NULL){
        aStr.push_back(pch);
        pch = strtok(NULL," ,");
      }
    }
//    for(unsigned int istr=0;istr<aStr.size();istr++){
//      std::cout << aStr[istr] << std::endl;
//    }
    
    std::vector< std::pair<Cad::CURVE_TYPE, std::vector<double> > > aVal;
    unsigned int istr  = 0;
    Com::CVector2D cur(0,0);
    for(;;){
      assert( !aStr[istr].empty() );      
      assert( aStr[istr][0] >= 0x3a );
      char ctl = aStr[istr][0];
      if( ctl == 'M' || ctl == 'm' ){
        bool is_rel = (ctl == 'm');        
        assert( aVal.size() == 0 );
        {
          std::vector<double> tmp;
          tmp.push_back( atof(aStr[istr+1].c_str()) );
          tmp.push_back( atof(aStr[istr+2].c_str()) );
          cur = Com::CVector2D(tmp[0],tmp[1]);
          aVal.push_back( std::make_pair(CURVE_END_POINT,tmp) );
          istr += 3;
        }
        for(;;){
          if( aStr[istr][0] >= 0x3a ) break;          
          std::vector<double> tmp;
          tmp.push_back( atof(aStr[istr+0].c_str()) );
          tmp.push_back( atof(aStr[istr+1].c_str()) );
          if( is_rel ){ tmp[0] += cur.x; tmp[1] += cur.y; }          
          aVal.push_back( std::make_pair(CURVE_LINE,tmp) ); 
          cur = Com::CVector2D(tmp[0],tmp[1]);
          istr += 2;
        }                
        continue;
      }
      else if( ctl == 'C' || ctl == 'c' ){
        bool is_rel = (ctl == 'c');
        istr += 1;
        for(;;){
          std::vector<double> tmp;
          for(unsigned int i=0;i<6;i++){
            tmp.push_back( atof(aStr[istr+i].c_str()) );
          }          
          if( is_rel ){ 
            tmp[0] += cur.x; tmp[1] += cur.y; 
            tmp[2] += cur.x; tmp[3] += cur.y; 
            tmp[4] += cur.x; tmp[5] += cur.y;             
          }
          aVal.push_back( std::make_pair(CURVE_BEZIER,tmp) );          
          cur = Com::CVector2D(tmp[4],tmp[5]);          
          istr += 6;
          if( aStr[istr][0] >= 0x3a ){
            break;          
          }
        }
        continue;
      }
      else if( ctl == 'L' || ctl == 'l' ){
        bool is_rel = (ctl == 'l');
        istr += 1;
        for(;;){
          std::vector<double> tmp;
          tmp.push_back( atof(aStr[istr+0].c_str()) );
          tmp.push_back( atof(aStr[istr+1].c_str()) );
          if( is_rel ){ tmp[0] += cur.x; tmp[1] += cur.y; }
          aVal.push_back( std::make_pair(CURVE_LINE,tmp) );          
          cur = Com::CVector2D(tmp[0],tmp[1]);          
          istr += 2;
          if( aStr[istr][0] >= 0x3a ){
            break;
          }
        }        
      }
      else if( ctl == 'z' ){
        assert( istr == aStr.size()-1 );
        break;
      }      
    }
    unsigned int id_l = cad_2d.AddLoop(aVal,0,scale).id_l_add;
    std::cout << "added loop: " << id_l << std::endl;
    res.push_back(id_l);
    ////
  }
  return res;
}

