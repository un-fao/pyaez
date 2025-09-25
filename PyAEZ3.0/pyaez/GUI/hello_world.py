import streamlit as st

st.title("Hello World!")

st.write("Welcome to your first Streamlit app.")

if st.button("Quit"):
    st.warning("Streamlit apps can't be 'quit' like desktop apps, but you can close the browser tab.")
