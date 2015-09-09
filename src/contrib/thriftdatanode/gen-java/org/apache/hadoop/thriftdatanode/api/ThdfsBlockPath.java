/**
 * Autogenerated by Thrift Compiler (0.7.0)
 *
 * DO NOT EDIT UNLESS YOU ARE SURE THAT YOU KNOW WHAT YOU ARE DOING
 */
package org.apache.hadoop.thriftdatanode.api;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.EnumMap;
import java.util.Set;
import java.util.HashSet;
import java.util.EnumSet;
import java.util.Collections;
import java.util.BitSet;
import java.nio.ByteBuffer;
import java.util.Arrays;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ThdfsBlockPath implements org.apache.thrift.TBase<ThdfsBlockPath, ThdfsBlockPath._Fields>, java.io.Serializable, Cloneable {
  private static final org.apache.thrift.protocol.TStruct STRUCT_DESC = new org.apache.thrift.protocol.TStruct("ThdfsBlockPath");

  private static final org.apache.thrift.protocol.TField LOCAL_BLOCK_PATH_FIELD_DESC = new org.apache.thrift.protocol.TField("localBlockPath", org.apache.thrift.protocol.TType.STRING, (short)1);
  private static final org.apache.thrift.protocol.TField LOCAL_META_PATH_FIELD_DESC = new org.apache.thrift.protocol.TField("localMetaPath", org.apache.thrift.protocol.TType.STRING, (short)2);

  public String localBlockPath; // required
  public String localMetaPath; // required

  /** The set of fields this struct contains, along with convenience methods for finding and manipulating them. */
  public enum _Fields implements org.apache.thrift.TFieldIdEnum {
    LOCAL_BLOCK_PATH((short)1, "localBlockPath"),
    LOCAL_META_PATH((short)2, "localMetaPath");

    private static final Map<String, _Fields> byName = new HashMap<String, _Fields>();

    static {
      for (_Fields field : EnumSet.allOf(_Fields.class)) {
        byName.put(field.getFieldName(), field);
      }
    }

    /**
     * Find the _Fields constant that matches fieldId, or null if its not found.
     */
    public static _Fields findByThriftId(int fieldId) {
      switch(fieldId) {
        case 1: // LOCAL_BLOCK_PATH
          return LOCAL_BLOCK_PATH;
        case 2: // LOCAL_META_PATH
          return LOCAL_META_PATH;
        default:
          return null;
      }
    }

    /**
     * Find the _Fields constant that matches fieldId, throwing an exception
     * if it is not found.
     */
    public static _Fields findByThriftIdOrThrow(int fieldId) {
      _Fields fields = findByThriftId(fieldId);
      if (fields == null) throw new IllegalArgumentException("Field " + fieldId + " doesn't exist!");
      return fields;
    }

    /**
     * Find the _Fields constant that matches name, or null if its not found.
     */
    public static _Fields findByName(String name) {
      return byName.get(name);
    }

    private final short _thriftId;
    private final String _fieldName;

    _Fields(short thriftId, String fieldName) {
      _thriftId = thriftId;
      _fieldName = fieldName;
    }

    public short getThriftFieldId() {
      return _thriftId;
    }

    public String getFieldName() {
      return _fieldName;
    }
  }

  // isset id assignments

  public static final Map<_Fields, org.apache.thrift.meta_data.FieldMetaData> metaDataMap;
  static {
    Map<_Fields, org.apache.thrift.meta_data.FieldMetaData> tmpMap = new EnumMap<_Fields, org.apache.thrift.meta_data.FieldMetaData>(_Fields.class);
    tmpMap.put(_Fields.LOCAL_BLOCK_PATH, new org.apache.thrift.meta_data.FieldMetaData("localBlockPath", org.apache.thrift.TFieldRequirementType.DEFAULT, 
        new org.apache.thrift.meta_data.FieldValueMetaData(org.apache.thrift.protocol.TType.STRING)));
    tmpMap.put(_Fields.LOCAL_META_PATH, new org.apache.thrift.meta_data.FieldMetaData("localMetaPath", org.apache.thrift.TFieldRequirementType.DEFAULT, 
        new org.apache.thrift.meta_data.FieldValueMetaData(org.apache.thrift.protocol.TType.STRING)));
    metaDataMap = Collections.unmodifiableMap(tmpMap);
    org.apache.thrift.meta_data.FieldMetaData.addStructMetaDataMap(ThdfsBlockPath.class, metaDataMap);
  }

  public ThdfsBlockPath() {
  }

  public ThdfsBlockPath(
    String localBlockPath,
    String localMetaPath)
  {
    this();
    this.localBlockPath = localBlockPath;
    this.localMetaPath = localMetaPath;
  }

  /**
   * Performs a deep copy on <i>other</i>.
   */
  public ThdfsBlockPath(ThdfsBlockPath other) {
    if (other.isSetLocalBlockPath()) {
      this.localBlockPath = other.localBlockPath;
    }
    if (other.isSetLocalMetaPath()) {
      this.localMetaPath = other.localMetaPath;
    }
  }

  public ThdfsBlockPath deepCopy() {
    return new ThdfsBlockPath(this);
  }

  @Override
  public void clear() {
    this.localBlockPath = null;
    this.localMetaPath = null;
  }

  public String getLocalBlockPath() {
    return this.localBlockPath;
  }

  public ThdfsBlockPath setLocalBlockPath(String localBlockPath) {
    this.localBlockPath = localBlockPath;
    return this;
  }

  public void unsetLocalBlockPath() {
    this.localBlockPath = null;
  }

  /** Returns true if field localBlockPath is set (has been assigned a value) and false otherwise */
  public boolean isSetLocalBlockPath() {
    return this.localBlockPath != null;
  }

  public void setLocalBlockPathIsSet(boolean value) {
    if (!value) {
      this.localBlockPath = null;
    }
  }

  public String getLocalMetaPath() {
    return this.localMetaPath;
  }

  public ThdfsBlockPath setLocalMetaPath(String localMetaPath) {
    this.localMetaPath = localMetaPath;
    return this;
  }

  public void unsetLocalMetaPath() {
    this.localMetaPath = null;
  }

  /** Returns true if field localMetaPath is set (has been assigned a value) and false otherwise */
  public boolean isSetLocalMetaPath() {
    return this.localMetaPath != null;
  }

  public void setLocalMetaPathIsSet(boolean value) {
    if (!value) {
      this.localMetaPath = null;
    }
  }

  public void setFieldValue(_Fields field, Object value) {
    switch (field) {
    case LOCAL_BLOCK_PATH:
      if (value == null) {
        unsetLocalBlockPath();
      } else {
        setLocalBlockPath((String)value);
      }
      break;

    case LOCAL_META_PATH:
      if (value == null) {
        unsetLocalMetaPath();
      } else {
        setLocalMetaPath((String)value);
      }
      break;

    }
  }

  public Object getFieldValue(_Fields field) {
    switch (field) {
    case LOCAL_BLOCK_PATH:
      return getLocalBlockPath();

    case LOCAL_META_PATH:
      return getLocalMetaPath();

    }
    throw new IllegalStateException();
  }

  /** Returns true if field corresponding to fieldID is set (has been assigned a value) and false otherwise */
  public boolean isSet(_Fields field) {
    if (field == null) {
      throw new IllegalArgumentException();
    }

    switch (field) {
    case LOCAL_BLOCK_PATH:
      return isSetLocalBlockPath();
    case LOCAL_META_PATH:
      return isSetLocalMetaPath();
    }
    throw new IllegalStateException();
  }

  @Override
  public boolean equals(Object that) {
    if (that == null)
      return false;
    if (that instanceof ThdfsBlockPath)
      return this.equals((ThdfsBlockPath)that);
    return false;
  }

  public boolean equals(ThdfsBlockPath that) {
    if (that == null)
      return false;

    boolean this_present_localBlockPath = true && this.isSetLocalBlockPath();
    boolean that_present_localBlockPath = true && that.isSetLocalBlockPath();
    if (this_present_localBlockPath || that_present_localBlockPath) {
      if (!(this_present_localBlockPath && that_present_localBlockPath))
        return false;
      if (!this.localBlockPath.equals(that.localBlockPath))
        return false;
    }

    boolean this_present_localMetaPath = true && this.isSetLocalMetaPath();
    boolean that_present_localMetaPath = true && that.isSetLocalMetaPath();
    if (this_present_localMetaPath || that_present_localMetaPath) {
      if (!(this_present_localMetaPath && that_present_localMetaPath))
        return false;
      if (!this.localMetaPath.equals(that.localMetaPath))
        return false;
    }

    return true;
  }

  @Override
  public int hashCode() {
    return 0;
  }

  public int compareTo(ThdfsBlockPath other) {
    if (!getClass().equals(other.getClass())) {
      return getClass().getName().compareTo(other.getClass().getName());
    }

    int lastComparison = 0;
    ThdfsBlockPath typedOther = (ThdfsBlockPath)other;

    lastComparison = Boolean.valueOf(isSetLocalBlockPath()).compareTo(typedOther.isSetLocalBlockPath());
    if (lastComparison != 0) {
      return lastComparison;
    }
    if (isSetLocalBlockPath()) {
      lastComparison = org.apache.thrift.TBaseHelper.compareTo(this.localBlockPath, typedOther.localBlockPath);
      if (lastComparison != 0) {
        return lastComparison;
      }
    }
    lastComparison = Boolean.valueOf(isSetLocalMetaPath()).compareTo(typedOther.isSetLocalMetaPath());
    if (lastComparison != 0) {
      return lastComparison;
    }
    if (isSetLocalMetaPath()) {
      lastComparison = org.apache.thrift.TBaseHelper.compareTo(this.localMetaPath, typedOther.localMetaPath);
      if (lastComparison != 0) {
        return lastComparison;
      }
    }
    return 0;
  }

  public _Fields fieldForId(int fieldId) {
    return _Fields.findByThriftId(fieldId);
  }

  public void read(org.apache.thrift.protocol.TProtocol iprot) throws org.apache.thrift.TException {
    org.apache.thrift.protocol.TField field;
    iprot.readStructBegin();
    while (true)
    {
      field = iprot.readFieldBegin();
      if (field.type == org.apache.thrift.protocol.TType.STOP) { 
        break;
      }
      switch (field.id) {
        case 1: // LOCAL_BLOCK_PATH
          if (field.type == org.apache.thrift.protocol.TType.STRING) {
            this.localBlockPath = iprot.readString();
          } else { 
            org.apache.thrift.protocol.TProtocolUtil.skip(iprot, field.type);
          }
          break;
        case 2: // LOCAL_META_PATH
          if (field.type == org.apache.thrift.protocol.TType.STRING) {
            this.localMetaPath = iprot.readString();
          } else { 
            org.apache.thrift.protocol.TProtocolUtil.skip(iprot, field.type);
          }
          break;
        default:
          org.apache.thrift.protocol.TProtocolUtil.skip(iprot, field.type);
      }
      iprot.readFieldEnd();
    }
    iprot.readStructEnd();

    // check for required fields of primitive type, which can't be checked in the validate method
    validate();
  }

  public void write(org.apache.thrift.protocol.TProtocol oprot) throws org.apache.thrift.TException {
    validate();

    oprot.writeStructBegin(STRUCT_DESC);
    if (this.localBlockPath != null) {
      oprot.writeFieldBegin(LOCAL_BLOCK_PATH_FIELD_DESC);
      oprot.writeString(this.localBlockPath);
      oprot.writeFieldEnd();
    }
    if (this.localMetaPath != null) {
      oprot.writeFieldBegin(LOCAL_META_PATH_FIELD_DESC);
      oprot.writeString(this.localMetaPath);
      oprot.writeFieldEnd();
    }
    oprot.writeFieldStop();
    oprot.writeStructEnd();
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder("ThdfsBlockPath(");
    boolean first = true;

    sb.append("localBlockPath:");
    if (this.localBlockPath == null) {
      sb.append("null");
    } else {
      sb.append(this.localBlockPath);
    }
    first = false;
    if (!first) sb.append(", ");
    sb.append("localMetaPath:");
    if (this.localMetaPath == null) {
      sb.append("null");
    } else {
      sb.append(this.localMetaPath);
    }
    first = false;
    sb.append(")");
    return sb.toString();
  }

  public void validate() throws org.apache.thrift.TException {
    // check for required fields
  }

  private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
    try {
      write(new org.apache.thrift.protocol.TCompactProtocol(new org.apache.thrift.transport.TIOStreamTransport(out)));
    } catch (org.apache.thrift.TException te) {
      throw new java.io.IOException(te);
    }
  }

  private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException {
    try {
      read(new org.apache.thrift.protocol.TCompactProtocol(new org.apache.thrift.transport.TIOStreamTransport(in)));
    } catch (org.apache.thrift.TException te) {
      throw new java.io.IOException(te);
    }
  }

}

